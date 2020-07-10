import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

def tagtracker(s, libhits, graphML, source, tag):
    
    for v in graphML.nodes():
        for j in libhits.columns:
            if int(v) in libhits['#Scan#'].tolist():
                graphML.nodes[v][j] = str(libhits[libhits['#Scan#'] == int(v)][j].iloc[-1])
            else:
                graphML.nodes[v][j] = 'nan'
            
    sel = s[s[tag].str.match(source, na = False)]['InChI_Key_Planar'].tolist()
    
    # select all nodes with library matches within selecte source category
    sel_ids = list()
    
    for v in graphML.nodes():
        if graphML.nodes[v]['InChIKey-Planar'] in sel:
            graphML.nodes[v]['TAGS_source'] = str(s[s['InChI_Key_Planar'] == graphML.nodes[v]['InChIKey-Planar']][tag].iloc[-1])
            #print(graphML.nodes[v]['TAGS_source'])
            tags = ['tag_' + s for s in graphML.nodes[v]['TAGS_source'].split('|')]
            for t in tags:
                graphML.nodes[v][t] = 1.0
            sel_ids.append(v)

    if sel_ids == []:
        print('No library hits contained the tags you are looking for')
    
    else: 
        # select all neighbours of the library matches
        sel_ids_neigh = list()
        
        for t in sel_ids:
            sel_ids_neigh.append([n for n in graphML.neighbors(str(t))])
            
        sel_ids_neigh = list(set([item for sublist in sel_ids_neigh for item in sublist]))
        
        # combine nodes with library matches and neighbours and create subgraph
        sel_ids_final = list(set(sel_ids_neigh + sel_ids))
                
        H = graphML.subgraph(sel_ids_final)
        
        # Plot Delta MZ histogram
        mdd = [item[2] for item in list(H.edges(data=True))]
        md = [d["mass_difference"] for d in mdd if "mass_difference" in d]
        md = [float(i) for i in md]
        md = [abs(round(elem, 2)) for elem in md]
        
        n, bins, patches = plt.hist(x=np.asarray(md), bins=100, color='#67a9cf',
                                alpha=0.7, rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('MZ Delta')
        plt.ylabel('Pairs')
        plt.title('Subnetwork MZ Delta Histogram')
    
        return H


def sourcetracker(s,libhits,column):
    
    tags = list()
    for v in libhits['InChIKey-Planar']:
        if v in s['InChI_Key_Planar'].tolist():
            tags.append(s[s['InChI_Key_Planar']== v][column].iloc[-1])
    
    tags_single = [word for line in tags for word in line.split('|')]
    
    # select metadata column to be plotted
    recounted = Counter(tuple(tags_single))

    # Make fake dataset
    height = [ v for v in recounted.values() ]
    bars = tuple([k for k in recounted])
    y_pos = np.arange(len(bars))

    # Create horizontal bars
    plt.barh(y_pos, height)

    # Create names on the y-axis
    plt.yticks(y_pos, bars)

    plt.xlabel('Number of Annotations')
    plt.ylabel(column)

    # Show graphic
    plt.show()
