#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
#%load_ext rpy2.ipython

sample_one = anndata.read_loom("/media/user/sdh/cmq_velo/out_220718/SRR_ALL.loom")


# In[2]:



sample_obs = pd.read_csv("/media/user/sdh/cmq_velo/out_220718/cellID_obs_T2D.csv")
umap_cord = pd.read_csv("/media/user/sdh/cmq_velo/out_220718/cell_embeddings_T2D.csv")
cell_clusters = pd.read_csv("/media/user/sdh/cmq_velo/out_220718/clusters_T2D.csv")
sample_one.uns['cell_clusters_colors']=np.array(["#8C564BFF" ,"#E377C2FF", "#C49C94FF", "#BCBD22FF", "#17BECFFF",
                 "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF" ,"#D62728FF", "#9467BDFF"])


# In[3]:


sample_one


# In[4]:


sample_obs


# In[5]:


sample_one.obs


# In[6]:


sample_obs


# In[7]:


#np.isin(a,b) 用于判定a中的元素在b中是否出现过，如果出现过返回True,否则返回False,最终结果为一个形状和a一模一样的数组。
sample_one = sample_one[np.isin(sample_one.obs.index,sample_obs["Cells(endo_as_sce)"])]


# In[8]:


sample_one.obs.index


# In[9]:


sample_one.obs.index


# In[10]:


umap = pd.read_csv("/media/user/sdh/cmq_velo/out_220718/cell_embeddings_T2D.csv")
umap


# In[11]:


sample_one_index = pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {0:'CellID'})


# In[12]:


sample_one_index


# In[13]:


umap = umap.rename(columns = {'Unnamed: 0':'CellID'})
cell_clusters =cell_clusters.rename(columns = {'Unnamed: 0':'CellID'})


# In[14]:


cell_clusters


# In[15]:


umap


# In[16]:


umap_ordered = sample_one_index.merge(umap,on='CellID')


# In[17]:


umap_ordered


# In[18]:


umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
sample_one.obsm['X_umap']


# In[19]:


cell_clusters_orderd =sample_one_index.merge(cell_clusters,on='CellID')


# In[20]:


cell_clusters


# In[21]:


cell_clusters_orderd
cell_clusters_orderd_1 = cell_clusters_orderd.iloc[:,1:]
cell_clusters_orderd_1.values
sample_one.obs['cell_clusters']= cell_clusters_orderd_1.values


# In[22]:


cell_clusters_orderd


# In[23]:


scv.pp.filter_and_normalize(sample_one,min_shared_counts=20, n_top_genes=3000)
#scv.pp.highly_variable_genes(sample_one)
scv.pp.moments(sample_one)
scv.tl.velocity(sample_one, mode = "stochastic")
scv.tl.velocity_graph(sample_one)
scv.pl.velocity_embedding(sample_one, basis='umap',color = ['cell_clusters'])


# In[24]:


scv.tl.velocity(sample_one)


# In[25]:


scv.tl.velocity_graph(sample_one)


# In[26]:


scv.pl.velocity_embedding_stream(sample_one, basis='umap',color = ['cell_clusters'])


# In[36]:


get_ipython().magic(u'pinfo scv.pl.velocity_embedding_stream')


# In[28]:


import os
print(os.path.abspath('.'))


# In[49]:


scv.pl.velocity_embedding_stream(sample_one, basis='umap',figsize=(10,10),color = ['cell_clusters'], 
                                 save='scv.pl.velocity_embedding_stream_T2D.svg')


# In[30]:


scv.pl.velocity_embedding_stream(sample_one, basis='umap',figsize=(10,10),
                                 color = ['cell_clusters'],arrow_size=3,arrow_style='->',linewidth=2, 
                                 save='scv.pl.velocity_embedding_stream_T2D_2.pdf' )
scv.pl.velocity_embedding_stream(sample_one, basis='umap',figsize=(10,10),
                                 color = ['cell_clusters'],arrow_size=2,linewidth=2, 
                                 save='scv.pl.velocity_embedding_stream_T2D_3.pdf' )


# In[80]:


scv.pl.velocity_embedding_grid(sample_one,layer=['velocity', 'spliced'],alpha =0.08,
                               figsize=(5.5,5), arrow_size=3,color = ['cell_clusters'],
                               save='scv.pl.velocity_embedding_grid_T2D_0805.pdf')


# In[57]:





# In[ ]:


scv.tl.rank_velocity_genes(sample_one, groupby='cell_clusters', min_corr=.3)

df = scv.DataFrame(sample_one.uns['rank_velocity_genes']['names'])
df.to_csv('rank_gene_ND.csv')
df.head()


# In[27]:


scv.tl.paga(sample_one, groups='cell_clusters')
df = scv.get_df(sample_one, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[76]:


scv.pl.paga(sample_one, basis='umap', size=660, alpha=.0,
            min_edge_width=2, node_size_scale=5,figsize=(5,5),save='scv.pl.paga_T2D_5.pdf')


# In[38]:


sample_obs_no4 = pd.read_csv("/media/user/sdh/cmq_velo/out_220718/cellID_obs_T2D_NO4.csv")
sample_one_no4= sample_one[np.isin(sample_one.obs.index,sample_obs_no4["Cells(endo_as_sce)"])]


# In[43]:


sample_one.uns['paga']


# In[39]:


sample_one_no4


# In[36]:


scv.pl.paga(sample_one_nd, basis='umap', size=30, alpha=.5,
            min_edge_width=2, node_size_scale=1.5,figsize=(6,6))


# In[45]:


scv.pl.proportions(sample_one,groupby= 'cell_clusters')


# In[46]:


scv.pl.velocity(sample_one, ['INS','CHL1','SLC2A2','MAFA','RFX6','PDX1'],figsize=(6,7),
                ncols=1,color = ['cell_clusters'],save='velocity_gene_T2D.pdf')


# In[47]:


scv.tl.rank_velocity_genes(sample_one, groupby='cell_clusters', min_corr=.3)

df = scv.DataFrame(sample_one.uns['rank_velocity_genes']['names'])
df.head()


# In[48]:


scv.tl.score_genes_cell_cycle(sample_one)
scv.pl.scatter(sample_one, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])


# In[49]:


scv.pl.velocity_graph(sample_one, threshold=.1)


# In[50]:


scv.tl.velocity_pseudotime(sample_one)
scv.pl.scatter(sample_one, color='velocity_pseudotime', cmap='gnuplot')


# In[51]:


scv.tl.recover_dynamics(sample_one)

scv.tl.latent_time(sample_one)
scv.pl.heatmap(sample_one, var_names=['INS','CHL1','SLC2A2','MAFA','RFX6','PDX1'], sortby='latent_time',
               col_color='cell_clusters', n_convolve=100)

