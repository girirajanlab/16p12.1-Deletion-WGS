import pandas as pd
import numpy as np

import networkx as nx

from itertools import permutations

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle, Wedge, Patch
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.colors as mcolors

matplotlib.rcParams['pdf.fonttype'] = 42

# Identify all parents in a generation
def find_parents(ped, mother_col='Mother', father_col='Father'):
	parent_pairs=ped[[mother_col, father_col]].drop_duplicates()
	parent_pairs=parent_pairs[~((parent_pairs[[mother_col, father_col]]=='0').all(axis=1))]

	# Return lists of tuples of parent pairs
	par_tup=[tuple(row) for row in parent_pairs[[mother_col, father_col]].to_numpy()]
	return(par_tup)

# Identify siblings
def find_siblings(ped, sample, half=False, mother_col='Mother', father_col='Father', iid_col='SGCode'):
	mother=ped[ped[iid_col]==sample][mother_col].to_list()[0]
	father=ped[ped[iid_col]==sample][father_col].to_list()[0]
	
	sibs=ped[(ped[iid_col]!=sample) & (ped[mother_col]==mother) & (ped[father_col]==father)][iid_col].to_list()
	
	if half:
		sibs+=ped[(ped[iid_col]!=sample) & ((ped[mother_col]==mother) | (ped[father_col]==father))][iid_col].to_list()
	
	return sibs

# Define the "generation" for each sample based on their parents/siblings/spouses/partners
def length_longest_path_to_node(G, node):
    all_ancestors = nx.ancestors(G, node)
    all_ancestors.add(node)
    G_sub = G.subgraph(all_ancestors)
    return nx.dag_longest_path_length(G_sub)

# Annotate generations in family
def generations(ped, mother_col='Mother', father_col='Father', iid_col='SGCode'):
	net=nx.DiGraph()
	net.add_nodes_from(ped[iid_col].to_list())
	for col in [mother_col, father_col]:
		net.add_edges_from([tuple(row) for row in ped[~ped[col].str.contains('UNRECRUITED')][[col, iid_col]].to_numpy()])
	gen_map={}
	for n in net.nodes():
		gen_map[n]=length_longest_path_to_node(net, n)
	
	# Ensure samples match their partners (use highest gen number)
	for n in net.nodes():
		pped=ped[(ped[[mother_col, father_col]]==n).any(axis=1)]
		if pped.shape[0]==0:
			continue
		max_gen=max([gen_map[i] for i in pped[mother_col].to_list()+pped[father_col].to_list() if 'UNRECRUITED' not in i])
		gen_map[n]=max_gen
	
	# Ensure each sample is one less than their parents
	for idx, row in ped.iterrows():
		n=row[iid_col]
		m=row[mother_col]
		f=row[father_col]
		parents=[gen_map[i] for i in [m,f] if 'UNRECRUITED' not in i]
		if len(parents)>0:
			gen_map[n]=min(parents)+1
	
	return ped[iid_col].map(gen_map)

# Annotate y position in pedigree based on generation
def y_position(generations, gen_spacing=0.2):
	return (-1)*generations*gen_spacing

# Function to score an order for TSP-like ordering
def score_order(order, reldf):
	dist_score=0
	n_score=0
	for i in range(len(order)):
		if i==len(order)-1:
			continue
		o1=order[i]
		o2=order[i+1]
		
		dist,n=list(reldf[(reldf[['S1', 'S2']].isin([o1, o2])).all(axis=1)][['shortest_path', 'n_shortest_paths']].itertuples(index=False, name=None))[0]
		dist_score+=dist
		n_score+=n
	
	return (dist_score, n_score)

# Define brute firce method for TSP-like ordering
def brute_tsp(net, samples):
	relatedness_lst=[]
	for s1 in samples:
		idx1=samples.index(s1)
		for s2 in samples:
			idx2=samples.index(s2)
			if idx1>=idx2:
				continue
			shortest_length=nx.shortest_path_length(net, s1, s2)
			num_shortest=len(list(nx.all_shortest_paths(net, s1, s2)))
			relatedness_lst.append([s1, s2, shortest_length, num_shortest])
	reldf=pd.DataFrame(relatedness_lst, columns=['S1', 'S2', 'shortest_path', 'n_shortest_paths'])
	
	# Create all possible orders
	orders=permutations(samples)
	ordf=pd.DataFrame({'Order':orders})
	ordf['dist_score'], ordf['n_score']=zip(*ordf.Order.apply(lambda x: score_order(x, reldf)))
	
	# Sort to find best order
	ordf['S1']=ordf.Order.str[0]
	ordf.sort_values(by=['dist_score', 'n_score', 'S1'], ascending=[True, False, True], inplace=True)
	order=ordf.Order.to_list()[0]
	
	return list(order)

# Define x positions
def x_position(ped, mother_col='Mother', father_col='Father', iid_col='SGCode', skip_symbol='0', partner_width=0.2, sib_width=0.2, unrel_width=0.4):
	# Create a network for the family
	net=nx.Graph()
	net.add_nodes_from(ped[iid_col].to_list())
	for col in [mother_col, father_col]:
		net.add_edges_from([tuple(row) for row in ped[~(ped[col].str.contains(skip_symbol))][[col, iid_col]].to_numpy()])
	
	# Define x positions by generation
	pos_map={}
	max_gen=ped.generations.max()
	for i in range(max_gen, -1, -1):
		# If this is the first generation, order individuals by relatedness
		if i==max_gen:
			gen_samps=ped[ped.generations==i][iid_col].to_list()
			order=brute_tsp(net, gen_samps)
			
			# Assign x values starting with first in order and adding sib_width if pair shares a parent and unrel_width if they do not
			x=0
			for i in range(len(order)):
				o=order[i]
				pos_map[o]=x
				
				if i==(len(order)-1):
					continue
				nexto=order[i+1]
				oped=ped[ped[iid_col].isin([o, nexto])]
				if len(list(oped[mother_col].unique()))<2 or len(list(oped[father_col].unique()))<2:
					x+=sib_width
				else:
					x+=unrel_width
		
		# For all other generations, base locations on where children are
		# If a sample has no children, base location on where siblings are
		else:
			gendf=ped[ped.generations==i].copy()
			if gendf.shape[0]==0:
				continue
			gendf['parent']=gendf[iid_col].isin(ped[mother_col].to_list()+ped[father_col].to_list())
			
			# Identify parent pairs
			par_tups=find_parents(ped[(ped[[mother_col, father_col]].isin(gendf[gendf.parent][iid_col].to_list()).any(axis=1))])
			
			# Check for samples in multiple pairs
			gendf['multi_pair']=gendf[iid_col].apply(lambda x: len([i for i in par_tups if x in i and '0' not in i])>1)
			
			# For parents in multiple pairs, set parent at the midpoint of all their children
			mpp_tups=[]
			mp_pars=gendf[gendf.multi_pair][iid_col].to_list()
			for par in mp_pars:
				kids=ped[(ped[[mother_col, father_col]]==par).any(axis=1)][iid_col].to_list()
				kid_pos=[pos_map[i] for i in kids]
				pos_map[par]=np.mean(kid_pos)
				
				# Make a list of parent tuples with these partners for later
				mpp_tups+=[i for i in par_tups if par in i]
				
			# For partners of parents with multiple pairs, set them partner_width away closer to their children
			for mpp in mpp_tups:
				if mpp[0] in mp_pars and mpp[1] in mp_pars:
					continue
				par=[i for i in mpp if i not in mp_pars][0]
				mp_par=[i for i in mpp if i in mp_pars][0]
				kids=ped[(ped[[mother_col, father_col]]==par).any(axis=1)][iid_col].to_list()
				kid_pos=[pos_map[i] for i in kids]
				mp_x=pos_map[mp_par]
				
				dir=1
				if np.mean(kid_pos)<mp_x:
					dir=-1
				
				pos_map[par]=mp_x+(dir*partner_width)
			
			# For partners in one pair, set them partner_width/2 away from the child mean and closer to siblings, if they have any
			for p in par_tups:
				if p in mpp_tups:
					continue
				# Check if they are part of a multiple partner pair (they may be included here if they are part of a pair with a missing sample)
				p1=p[0]
				p2=p[1]
				escape=False
				for mpp in mpp_tups:
					if p1 in mpp or p2 in mpp:
						escape=True
				if escape:
					continue
				
				# Check if either partner has a sibling
				if skip_symbol in p[0] or skip_symbol in p[1]:
					p1=[i for i in p if skip_symbol not in i][0]
					p2=[i for i in p if skip_symbol in i][0]
					par_order=[p1, p2]
				elif len(find_siblings(ped, p1))>0:
					par_order=[p1, p2]
				elif len(find_siblings(ped, p2))>0:
					par_order=[p2, p1]
				else:
					par_order=[p1, p2]

				dir=1
				for par in par_order:
					if skip_symbol in par:
						continue
					kids=ped[(ped[[mother_col, father_col]]==par).any(axis=1)][iid_col].to_list()
					kid_pos=[pos_map[i] for i in kids]
					kid_mean=np.mean(kid_pos)
					
					sibs=find_siblings(ped, par)
					if sibs!=[] and par==par_order[0]:
						sib_pos=[pos_map[i] for i in sibs if i in pos_map.keys()]
						if np.mean(sib_pos)<kid_mean:
							dir=-1
					
					pos_map[par]=kid_mean+(dir*(partner_width/2))
					dir=dir*-1
					# If a single parent, check if too close to another sample and flip as needed
					if skip_symbol in par_order[0] or skip_symbol in par_order[1]:
						taken_pos=[pos_map[i] for i in gendf[iid_col].to_list() if i!=par and i in pos_map.keys()]
						for tp in taken_pos:
							if round(tp, 3)>round(pos_map[par]-sib_width, 3) and round(tp, 3)<round(pos_map[par]+sib_width, 3):
								pos_map[par]=kid_mean+(dir*(partner_width/2))
								break
			
			# For individuals without children, set them sib_width away, shifting as needed to not overlap another sample
			siblings=[i for i in gendf[iid_col].to_list() if i not in pos_map.keys()]
			taken_pos=[pos_map[i] for i in gendf[iid_col].to_list() if i not in siblings]
			for sib in siblings:
				sibs=find_siblings(ped, sib)
				if sibs==[]:
					sibs=find_siblings(ped, sib, half=True)
				sib_pos=[pos_map[i] for i in sibs if i in pos_map.keys()]
				nonsib_pos=[pos_map[i] for i in pos_map.keys() if i not in sibs and i in gendf[iid_col].to_list()]
				if min(sib_pos)<min(nonsib_pos):
					dir=-1
					startx=max(sib_pos)
				else:
					dir=1
					startx=min(sib_pos)
				
				while True:
					# Make sure sibling is at least sib_width away from all other points
					minr=startx-(0.9*sib_width)
					maxr=startx+(0.9*sib_width)
					
					in_range=[i for i in taken_pos if i>=minr and i<=maxr]
					if len(in_range)==0:
						pos_map[sib]=startx
						taken_pos.append(startx)
						break
					else:
						startx+=(dir*sib_width)
				
	return ped[iid_col].map(pos_map)

# Add a patch representing a person
def plot_shape(ax, sex, x, y, shape_width=0.1, fc='None', ec='k', lw=1, zo=2):
	# location represents the center of the shape
	if sex=='M':
		ax.add_artist(Rectangle((x-(shape_width/2), y-(shape_width/2)), shape_width, shape_width, fc=fc, edgecolor=ec, lw=lw, zorder=zo))
	elif sex=='F':
		ax.add_artist(Circle((x, y), shape_width/2, fc=fc, edgecolor=ec, lw=lw, zorder=zo))
	else:
		ax.add_artist(Rectangle((x-(shape_width/2), y-(shape_width/2)), shape_width, shape_width, fc=fc, edgecolor=ec, lw=lw, zorder=zo, rotation_point='center', angle=45))

# Draw lines up from individuals to connect to parents
def add_lines_up(ax, x, y, shape_width=0.1, gen_spacing=0.2):
	ax.plot([x, x], [y+shape_width/2, y+gen_spacing/2], color='k', zorder=-1)

# Draw lines connecting siblings and parents to children
def add_connecting_lines(ax, ped, mother_col='Mother', father_col='Father', iid_col='SGCode', skip_symbol='0', sib_width=0.2, shape_width=0.1, gen_spacing=0.2, zo=-2):
	# Draw lines connecting family members
	gens=sorted(list(ped.generations.unique()))
	gens.reverse()
	for g in gens:
		subdf=ped[ped.generations==g].copy()
		# Connect children
		if g!=gens[-1]:
			# Find all parent pairs in above generation
			par_tup=find_parents(subdf)
			for pt in par_tup:
				sibdf=subdf[(subdf[[mother_col, father_col]].isin(pt)).all(axis=1)]
				if skip_symbol in pt[0] or skip_symbol in pt[1]:
					sibdf=sibdf[(~sibdf[mother_col].str.contains(skip_symbol)) | (~sibdf[father_col].str.contains(skip_symbol))]
				# Add a line connecting siblings
				if sibdf.shape[0]>1:
					ax.plot([sibdf.x.min(), sibdf.x.max()], [sibdf.y.min()+gen_spacing/2, sibdf.y.min()+gen_spacing/2], color='k', zorder=zo)
				# Add line from sibling middle to parent mid
				par_mid=ped[ped[iid_col].isin(pt)]['x'].mean()
				# If there is only one parent, draw a line to the sibling mid
				if skip_symbol in pt[0] or skip_symbol in pt[1]:
					par_mid=sibdf.x.mean()
				ax.plot([par_mid, par_mid], [sibdf.y.min()+gen_spacing/2, sibdf.y.min()+gen_spacing], color='k', zorder=zo)
				# If the parent line and sibling line don't match, draw a connecting line
				if par_mid!=sibdf.x.mean():
					ax.plot([par_mid, sibdf.x.mean()], [sibdf.y.min()+gen_spacing/2, sibdf.y.min()+gen_spacing/2], color='k', zorder=zo)
		# Connect partners
		if g!=gens[0]:
			# Get all partner pairs in this generation
			par_tup=find_parents(ped[ped.generations==(g+1)])
			for pt in par_tup:
				ptdf=subdf[(subdf[iid_col].isin(pt)) & (~subdf[iid_col].str.contains(skip_symbol))]
				# Draw line connecting partners
				if skip_symbol not in pt[0] and skip_symbol not in pt[1]:
					ax.plot([ptdf.x.min()+shape_width/2, ptdf.x.max()-shape_width/2], [ptdf.y.min(), ptdf.y.min()], color='k', zorder=zo)
				# If there is only one partner, connect line to child average
				else:
					kid_mean=ped[(ped.generations==(g+1)) & ((ped[[mother_col, father_col]].isin(pt)).all(axis=1))].x.mean()
					dir=1
					if kid_mean<ptdf.x.mean():
						dir=-1
					ax.plot([ptdf.x.min()+(dir*shape_width/2), kid_mean], [ptdf.y.min(), ptdf.y.min()], color='k', zorder=zo)

# Add color patches representing phenotypes
def plot_phenos(ax, sex, x, y, phenos, color_dict, shape_width=0.1, zo=0):
	# Add in blocks of color to represent phenotypes
	if sex=='M':
		# For males, draw rectangles with width=recnt_width/(number of phenotypes)
		w=shape_width/len(phenos)
		locx=x-(shape_width/2)
		for p in phenos:
			ax.add_artist(Rectangle((locx, y-(shape_width/2)), w, 0.1, facecolor=color_dict[p], edgecolor=None, zorder=zo))
			locx+=w
	elif sex=='F':
		# For females, draw pie slices for phenotypes
		d_theta=360/len(phenos)
		theta=90
		for p in phenos:
			ax.add_artist(Wedge((x, y), shape_width/2, theta, theta+d_theta, facecolor=color_dict[p], edgecolor=None, zorder=zo))
			theta+=d_theta

# Add a custom legend for phenotypes
def custom_pheno_legend(ax, color_dict):
	legend_elem=[]
	for p in color_dict.keys():
		legend_elem.append(Patch(fc=color_dict[p], edgecolor='None', label=p))
	legend_elem.append(Patch(fc='white', edgecolor='red', lw=2, label='16p12.1 deletion carrier'))
	ax.legend(handles=legend_elem, loc='center', ncols=2)

# Add sample IDs under shapes
def label_individuals(ax, ped, iid_col='SGCode', shape_width=0.1, gen_spacing=0.2, fs=7):
	for idx, row in ped.iterrows():
		if row[iid_col]!=row[iid_col]:
			continue
		ax.text(row.x, row.y-(shape_width/2)-(gen_spacing/5), row[iid_col], color='k', fontsize=fs, ha='center', va='bottom')

# Add custom legend for burden
def custom_burden_legend(ax, color_dict):
	legend_elem=[]
	for v in color_dict.keys():
		legend_elem.append(Line2D([0], [0], color=color_dict[v], label=v))
	ax.legend(handles=legend_elem, loc='center')

# placeholder text for spacing