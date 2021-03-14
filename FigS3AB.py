from rohan.global_imports import *

# read table
df2=read_table('dfrequency.pqt')
# filter table
df2=df2.rd.filter_rows({'mutation format':'amino acid'})
# plot 
def plot_reg_replicate(df):
    label=df.name
    dplot=df.pivot_table(index=['position reference','mutated'],columns='replicate#',values='frequency (log2 scale)')
    dplot.columns=[f"replicate{c}" for c in dplot]
    from rohan.dandage.plot.scatter import plot_reg
    from rohan.dandage.plot.ax_ import set_equallim
    xys=list(itertools.combinations(dplot.columns.tolist(),2))
    fig,axs=plt.subplots(nrows=1,ncols=len(xys),figsize=[3.5*len(xys),3])
    for ax,xy in zip(np.ravel(axs),xys):
        ax=plot_reg(dplot.reset_index(),xy[0],xy[1],
                scafmt='sca',trendline=True,
                params_scatter={'color':'k','alpha':0.5},
                ax=ax)
        ax=set_equallim(ax,diagonal=True)
        ax.grid(True)
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax.axis('scaled')
    savefig(f"plot/scatter_{' '.join(label)}",fmts=['svg'])
df2.loc[(df2['sample']=='ABP1') & (df2['condition'].isin(['HUA2_S2',
                                                          'LSB3_S2',
                                                          'HUA2_MTX2',
                                                          'LSB3_MTX2',
                                                         ])),:].groupby(['sample','condition']).apply(plot_reg_replicate)