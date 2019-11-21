import pickle
import networkx as nx
import re

G = None
with (open('non-closed_Directed_KG_triples_no-mesh.gpickle', 'rb')) as openfile:
    G = pickle.load(openfile)

# There are 109272 nodes
len(G.nodes)

# There are 1461016 edges
len(G.edges)

# When made undirected it is not fully connected (but pretty much is)
Gundir = G.to_undirected()
nx.is_connected(Gundir)

# There are 4 sub graphs, one of which has 109253 nodes...
gs = nx.connected_component_subgraphs(Gundir)
for i, sg in enumerate(gs):
    print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))

# Find the shortest path between two things
# HP:0007256
# SEMA4A gene
# http://purl.uniprot.org/geneid/64218
# The disease Rod-cone dystrophy Retinitis pigmentosa
# https://hpo.jax.org/app/browse/term/HP:0000510

#Use the gene as the source, pheno as the target
print(nx.shortest_path(G,source='<http://purl.uniprot.org/geneid/64218>', target='<http://purl.obolibrary.org/obo/HP_0000510>'))
# Find all Sema4a and RP paths
paths = nx.all_simple_paths(G, source='<http://purl.uniprot.org/geneid/64218>', target='<http://purl.obolibrary.org/obo/HP_0000510>')
#don't actually try to print all of them, there are way to many. I sat here for 10 minutes and never even got the lenght of the number of paths
# print(len(list(paths)))


print(nx.shortest_path(Gundir,target='<http://purl.obolibrary.org/obo/HP_0007256>', source='<https://reactome.org/content/detail/R-HSA-111469>'))
print(nx.shortest_path(G,target='<http://purl.obolibrary.org/obo/HP_0005978>', source='<https://reactome.org/content/detail/R-HSA-111469>'))


# Find all paths between two things
for path in nx.all_simple_paths(G, target='<http://purl.obolibrary.org/obo/DOID_9351>', source='<https://hpo.jax.org/app/browse/term/HP:0000007>'):
    print(path)


# There are 19 sources
# There are their url bases
# {'http://orcid.org/', 'http://purl.org/dc/elements/', 'https://github.com/obophenotype/uberon/wiki/', 'http://swrl.stanford.edu/ontologies/', 'http://www.ncbi.nlm.nih.gov/pubmed/', 'https://en.wikipedia.org/wiki/', 'http://www.w3.org/2000/01/', 'http://www.w3.org/2002/07/', 'http://www.geneontology.org/formats/', 'http://purl.obolibrary.org/obo/', 'http://purl.obolibrary.org/obo/ro/docs/', 'http://purl.obolibrary.org/obo/go/releases/', 'http://purl.obolibrary.org/obo/ro/', 'http://www.w3.org/2003/11/', 'http://purl.obolibrary.org/obo/hp/', 'http://ontologydesignpatterns.org/wiki/', 'https://reactome.org/content/detail/', 'http://purl.uniprot.org/geneid/', 'http://xmlns.com/foaf/'}
# This is the one wikipedia node:
# <https://en.wikipedia.org/wiki/Allen%27s_interval_algebra>
url_set = set()
for i in G.nodes:
    r1 = re.findall(r"https?://[a-zA-Z0-9.]*[\/a-zA-Z0-9_]*\/", i)
    if len(r1) > 0:
        if r1[0] =='https://en.wikipedia.org/wiki/':
            print(i)
        url_set.add(r1[0])

# Get the average out degree and in degree of all nodes
# Find the max
out_sum = 0
in_sum = 0
count = 0
best_in = 0
best_out = 0
best_in_id = ''
best_out_id = ''
for i in G.nodes:
    r1 = re.findall(r"w3\.org", i)
    #skip anything from w3 school, these are not medical nodes
    if len(r1) > 0:
        continue
    outD = G.out_degree(i)
    inD = G.in_degree(i)
    out_sum += outD
    in_sum += inD
    if inD > best_in:
        best_in = inD
        best_in_id = i
    if outD > best_out:
        best_out = outD
        best_out_id = i
    count += 1
avg_in = in_sum / count
avg_out = out_sum / count
print(avg_in)
# 12.682951830936931
print(avg_out)
# 13.371770348065642
print(best_in)
# 9407
print(best_in_id)
# protein binding
# <http://purl.obolibrary.org/obo/GO_0005515>
print(best_out)
# 1584
print(best_out_id)
# <http://purl.uniprot.org/geneid/7311>



Gundir2 = Gundir.copy()
# Remove nodes from w3.org from undirected graph
for i in Gundir.nodes:
    if len(re.findall(r"w3\.org", i)) > 0:
        Gundir2.remove_node(i)

# There are 2096 sub graphs, one of which has 106339 nodes...
# Most have just 1, 12 have more than 1
# subgraph 0 has 106339 nodes
# subgraph 3 has 505 nodes
# subgraph 11 has 295 nodes
# subgraph 83 has 3 nodes
# subgraph 119 has 11 nodes
# subgraph 136 has 8 nodes
# subgraph 192 has 2 nodes
# subgraph 321 has 4 nodes
# subgraph 480 has 2 nodes
# subgraph 553 has 3 nodes
# subgraph 746 has 2 nodes
# subgraph 1291 has 2 nodes
gs = nx.connected_component_subgraphs(Gundir2)
not_1 = 0
for i, sg in enumerate(gs):
    if sg.number_of_nodes() > 1:
        not_1 += 1
        print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
print(not_1)


# What is the compoisition of the 3 large graphs

def get_source_composition(g):
    m = {'http://orcid.org/':0, 'http://purl.org/dc/elements/':0, 'https://github.com/obophenotype/uberon/wiki/':0, 'http://swrl.stanford.edu/ontologies/':0, 'http://www.ncbi.nlm.nih.gov/pubmed/':0, 'https://en.wikipedia.org/wiki/':0, 'http://www.w3.org/2000/01/':0, 'http://www.w3.org/2002/07/':0, 'http://www.geneontology.org/formats/':0, 'http://purl.obolibrary.org/obo/':0, 'http://purl.obolibrary.org/obo/ro/docs/':0, 'http://purl.obolibrary.org/obo/go/releases/':0, 'http://purl.obolibrary.org/obo/ro/':0, 'http://www.w3.org/2003/11/':0, 'http://purl.obolibrary.org/obo/hp/':0, 'http://ontologydesignpatterns.org/wiki/':0, 'https://reactome.org/content/detail/':0, 'http://purl.uniprot.org/geneid/':0, 'http://xmlns.com/foaf/':0}
    for i in g.nodes:
        r1 = re.findall(r"https?://[a-zA-Z0-9.]*[\/a-zA-Z0-9_]*\/", i)
        if len(r1) > 0:
            m[r1[0]] += 1
    return m

# The subgraph with 505 and 295 nodes are made entirely of <'http://purl.obolibrary.org/obo/'>
# The large graphi is made mostly of <'http://purl.obolibrary.org/obo/'> with significant amounts of 'https://reactome.org/content/detail/' and 'http://purl.uniprot.org/geneid/'
# 0
# {'http://orcid.org/': 3, 'http://purl.org/dc/elements/': 0, 'https://github.com/obophenotype/uberon/wiki/': 1, 'http://swrl.stanford.edu/ontologies/': 0, 'http://www.ncbi.nlm.nih.gov/pubmed/': 5, 'https://en.wikipedia.org/wiki/': 1, 'http://www.w3.org/2000/01/': 0, 'http://www.w3.org/2002/07/': 0, 'http://www.geneontology.org/formats/': 1, 'http://purl.obolibrary.org/obo/': 75417, 'http://purl.obolibrary.org/obo/ro/docs/': 1, 'http://purl.obolibrary.org/obo/go/releases/': 1, 'http://purl.obolibrary.org/obo/ro/': 1, 'http://www.w3.org/2003/11/': 0, 'http://purl.obolibrary.org/obo/hp/': 0, 'http://ontologydesignpatterns.org/wiki/': 4, 'https://reactome.org/content/detail/': 11969, 'http://purl.uniprot.org/geneid/': 18935, 'http://xmlns.com/foaf/': 0}
# 3
# {'http://orcid.org/': 0, 'http://purl.org/dc/elements/': 0, 'https://github.com/obophenotype/uberon/wiki/': 0, 'http://swrl.stanford.edu/ontologies/': 0, 'http://www.ncbi.nlm.nih.gov/pubmed/': 0, 'https://en.wikipedia.org/wiki/': 0, 'http://www.w3.org/2000/01/': 0, 'http://www.w3.org/2002/07/': 0, 'http://www.geneontology.org/formats/': 0, 'http://purl.obolibrary.org/obo/': 505, 'http://purl.obolibrary.org/obo/ro/docs/': 0, 'http://purl.obolibrary.org/obo/go/releases/': 0, 'http://purl.obolibrary.org/obo/ro/': 0, 'http://www.w3.org/2003/11/': 0, 'http://purl.obolibrary.org/obo/hp/': 0, 'http://ontologydesignpatterns.org/wiki/': 0, 'https://reactome.org/content/detail/': 0, 'http://purl.uniprot.org/geneid/': 0, 'http://xmlns.com/foaf/': 0}
# 11
# {'http://orcid.org/': 0, 'http://purl.org/dc/elements/': 0, 'https://github.com/obophenotype/uberon/wiki/': 0, 'http://swrl.stanford.edu/ontologies/': 0, 'http://www.ncbi.nlm.nih.gov/pubmed/': 0, 'https://en.wikipedia.org/wiki/': 0, 'http://www.w3.org/2000/01/': 0, 'http://www.w3.org/2002/07/': 0, 'http://www.geneontology.org/formats/': 0, 'http://purl.obolibrary.org/obo/': 295, 'http://purl.obolibrary.org/obo/ro/docs/': 0, 'http://purl.obolibrary.org/obo/go/releases/': 0, 'http://purl.obolibrary.org/obo/ro/': 0, 'http://www.w3.org/2003/11/': 0, 'http://purl.obolibrary.org/obo/hp/': 0, 'http://ontologydesignpatterns.org/wiki/': 0, 'https://reactome.org/content/detail/': 0, 'http://purl.uniprot.org/geneid/': 0, 'http://xmlns.com/foaf/': 0}
gs = nx.connected_component_subgraphs(Gundir2)
for i, sg in enumerate(gs):
    if i == 11 or i == 3 or i ==0:
        print(i)
        print(get_source_composition(sg))

for i in Gundir.nodes:
    if len(re.findall(r"http://www\.geneontology\.org/formats/", i)) > 0:
        print(i)


http://purl.uniprot.org/geneid/64218