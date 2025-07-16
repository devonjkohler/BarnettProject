import networkx as nx

def create_graph_from_complexes(complex_list, complex_dictionary, n_tissues):
    """Create a NetworkX graph from a list of complexes and their gene members.
    
    Parameters
    ----------
    complex_list : list of str
        List of complex names to include in the graph.
    complex_dictionary : dict
        Dictionary mapping complex names to lists of gene names.
    n_tissues : int
        The number of tissues being analyzed, used to set initial edge weights.
        
    Returns
    -------
    G : networkx.Graph
        An undirected graph where nodes are genes and edges connect genes 
        that are part of the same complex.
    """

    subset_genes = {k: complex_dictionary[k] for k in complex_list}

    G = nx.Graph()

    for complex, genes in subset_genes.items():
        # Add edges between all pairs of genes in the complex
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                if genes[i] and genes[j]:  # Ensure neither gene is an empty string
                    G.add_edge(genes[i], genes[j], 
                               complex=complex, weight=n_tissues-1)
                    
    return G

def build_tissue_network(graph, da_results, tissue_label, pval_threshold=0.2):
    """Build a tissue-specific network from the given graph and differential abundance results.

    Parameters
    ----------
    graph : networkx.Graph
        The input graph representing Protein interactions.
    da_results : pd.DataFrame
        A DataFrame containing differential abundance results with columns:
        'Protein', 'pval', and 'logfc'.
    tissue_label : str
        The label of the tissue to build the network for (e.g., 'B', 'M').
    pval_threshold : float
        The p-value threshold for filtering differential abundance results.

    Returns
    -------
    networkx.Graph
        A subgraph of the input graph containing only the tissue-specific interactions.
    """
    tissue_graph = graph.copy()

    # Loop over nodes in G and remove all edges for nodes where bm_results p_value < 0.1
    for node in tissue_graph.nodes:
        node_results = da_results[da_results["Protein"] == node]
        
        for i in range(len(node_results)):
            comp = node_results["comparison"].values[i]
            pval = node_results["p_value"].values[i]
            
            # If the comparison is of the form "X-Y" and tissue_label is Y, flip logFC
            conds = comp.split("-")
            
            if conds[1] == tissue_label:
                # Flip logFC sign for tissue_label as the second condition
                logfc = -node_results["logFC"].values[i]
            else:
                logfc = node_results["logFC"].values[i]
            
            if (pval < pval_threshold) & (logfc < 0):
                # For all connecting edges, reduce weight by one
                for neighbor in tissue_graph.neighbors(node):
                    tissue_graph[node][neighbor]["weight"] -= 1
    
    edges_to_remove = [(u, v) for u, v, d in tissue_graph.edges(data=True) \
        if d.get("weight", 1) <= 0]
    tissue_graph.remove_edges_from(edges_to_remove)
  
    return tissue_graph