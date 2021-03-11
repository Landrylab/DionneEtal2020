import os
os.chdir("/Users/emilie/Dropbox/Maitrise/Microscopie")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import seaborn as sns
sns.set(style="white")


## LOAD AND MERGE DATA FROM ANALYZED FOLDERS

raw_trajectories = pd.DataFrame()
cellprops = pd.DataFrame()
FileNotFound = []

# make list of all the image IDs
dirs = [f for f in os.listdir('AnalysedImages') if not f.startswith(('.', 'cfg', 'segmentation_cell'))]

# for each image, append cell properties and GFP trajectories
for i in dirs:
    data_cellprops = [f for f in os.listdir(f'AnalysedImages/{i}') if f.endswith('.tif.segmented.tif.npy.cellprops.tsv')]
    data_cellprops = pd.read_csv(f'AnalysedImages/{i}/{data_cellprops[0]}', sep='\t')
    data_cellprops.rename(columns={'cell #':'cell'}, inplace=True)
    data_cellprops['prop'] = data_cellprops['image type'] + '-' + data_cellprops['property']
    data_cellprops = data_cellprops.pivot_table(index='cell', columns='prop', values='value')
    data_cellprops['imageID'] = i
    data_cellprops['strain'] = i.split('_')[0]
    cellprops = cellprops.append(data_cellprops)
    cells = [f for f in os.listdir(f'AnalysedImages/{i}/cells') if
             not f.startswith('.')]  # each cell has a separate folder
    for j in cells:
        try:
            data_dist = pd.read_csv(f'AnalysedImages/{i}/cells/{j}/d3filter_returns_distance.tsv', sep='\t')
            data_dist['imageID'] = i
            data_dist['strain'] = i.split('_')[0]
            data_dist['cell'] = int(j.replace('cell', ''))
            raw_trajectories = raw_trajectories.append(data_dist)
        except FileNotFoundError: # some cells have no dots detected
            FileNotFound.append(f'd3filter_returns_distance not found : {i} {j}')

# check the dataframes
raw_trajectories.info()
cellprops.info()

raw_trajectories.drop(columns=['Unnamed: 0', 'level_0', 'level_1'], inplace=True)
cellprops.reset_index(inplace=True)

# save raw dataframes
raw_trajectories.to_csv('SLA1GFPraw_trajectories.csv')
cellprops.to_csv('SLA1GFPcellprops.csv')


## FORMAT STRAIN NAMES

raw_names = []
new_names = []

for i in raw_trajectories['strain'].unique():
    raw_names.append(i)
    new = ''
    i = str(i).replace('WT', 'W')
    i = i.replace('x', 'D')
    for j in range(3):
        if i[j] == 'W':
            new += (str(j + 1))
        else:
            new += i[j]
    new = new[0] + '|' + new[1] + '|' + new[2]
    new_names.append(new)

# check the result
for i in range(len(raw_names)):
    print(f'{raw_names[i]} to {new_names[i]}')

# change names in dataframes
for i in [raw_trajectories, cellprops]:
    i.replace(raw_names, new_names, inplace=True)

# check the dataframes strain names
raw_trajectories['strain'].unique()
cellprops['strain'].unique()


## CONVERT UNITS

# distance units in nm
for i in ['distance delta', 'distance effective', 'distance effective from centroid',
          'distance effective from centroid max', 'distance effective from centroid min',
          'distance effective from centroid per frame', 'distance effective per frame',
          'distance total', 'distance total per frame']:
    raw_trajectories[i] = raw_trajectories[i]*43.33 # 43.33 nm/pixel

# add time and lifetime columns in seconds
raw_trajectories['time'] = raw_trajectories['frame'] - raw_trajectories['frame min'] # 1 frame/s
raw_trajectories['lifetime'] = raw_trajectories['frame max'] - raw_trajectories['frame min']

# check the updates
raw_trajectories.head()


## FIND AND FILTER MISDETECTED CELLS

# look at cells properties
# calculate difference between major and minor axis compared to the diameter
cellprops['diff_axis_length'] = cellprops['bright-major_axis_length'] - cellprops['bright-minor_axis_length']
cellprops['diff_axis/diameter'] = cellprops['diff_axis_length'] / cellprops['bright-equivalent_diameter']

# calculate circularity
cellprops['circularity'] = 4*np.pi*cellprops['bright-area'] / cellprops['bright-perimeter']**2

# calculate difference of mean intensity compared to the cells of the same image
cellprops.set_index('imageID', inplace=True)
cellprops['mean_image_intensity'] = 0
cellprops.update(cellprops.groupby(cellprops.index)['bright-mean_intensity'].mean().to_frame('mean_image_intensity'))
cellprops['diff_intensity'] = cellprops['bright-mean_intensity'] - cellprops['mean_image_intensity']
cellprops.reset_index(inplace=True)

# fix thresholds
area_min = cellprops['bright-area'].quantile(q=0.05) # small cells are frequently misdetected
diff_intensity_max = cellprops['diff_intensity'].quantile(q=0.85) # cells on top of each other are brighter
circularity_min = cellprops['circularity'].quantile(q=0.05) # deformed cells are misdetected
solidity_min = cellprops['bright-solidity'].quantile(q=0.05)  # deformed cells are misdetected
diff_axis_length_max = cellprops['diff_axis/diameter'].quantile(q=0.95)  # deformed cells are misdetected

# plot distributions with thresholds
fig, ax = plt.subplots(ncols=3, nrows=3, figsize=(15, 15))
plt.subplots_adjust(wspace = 0.3)

ax0 = sns.distplot(cellprops['bright-area'], ax=ax[0,0], color='black')
ax0.axvline(area_min, color='b')

ax1 = sns.distplot(cellprops['diff_intensity'], ax=ax[0,1], color='black')
ax1.axvline(diff_intensity_max, color='b')

ax2 = sns.scatterplot(x='bright-area', y='diff_intensity',
                      data=cellprops, marker='o', alpha=.2, ax=ax[0,2], color='black')
ax2.axvline(area_min, color='b')
ax2.axhline(diff_intensity_max, color='b')

ax3 = sns.distplot(cellprops['diff_axis/diameter'], ax=ax[1,0], color='black')
ax3.axvline(diff_axis_length_max, color='b')

ax4 = sns.distplot(cellprops['bright-solidity'], ax=ax[1,1], color='black')
ax4.axvline(solidity_min, color='b')

ax5 = sns.distplot(cellprops['circularity'], ax=ax[1,2], color='black')
ax5.axvline(circularity_min, color='b')

ax6 = sns.scatterplot(x='diff_axis/diameter', y='bright-solidity',
                      data=cellprops, marker='o', alpha=.2, ax=ax[2,0], color='black')
ax6.axvline(diff_axis_length_max, color='b')
ax6.axhline(solidity_min, color='b')

ax7 = sns.scatterplot(x='diff_axis/diameter', y='circularity',
                      data=cellprops, marker='o', alpha=.2, ax=ax[2,1], color='black')
ax7.axvline(diff_axis_length_max, color='b')
ax7.axhline(circularity_min, color='b')

ax8 = sns.scatterplot(x='bright-solidity', y='circularity',
                      data=cellprops, marker='o', alpha=.2, ax=ax[2,2], color='black')
ax8.axvline(solidity_min, color='b')
ax8.axhline(circularity_min, color='b')

plt.show()

# subset rejected cells
props_rejected = cellprops.loc[((cellprops['bright-area']<area_min)
                                | (cellprops['diff_intensity']>diff_intensity_max)
                                | (cellprops['circularity']<circularity_min)
                                | (cellprops['bright-solidity']<solidity_min)
                                | (cellprops['diff_axis/diameter']>diff_axis_length_max)), :]

# check number of cells
len(props_rejected) # 577 cells rejected / 2417

# confirm rejected cells by manually looking at some images
props_rejected.reset_index(inplace=True)
props_rejected.loc[(props_rejected['imageID'] == 'xxx_1202_002'), ['cell', 'bright-area', 'bright-mean_intensity',
                                                                   'circularity', 'bright-solidity',
                                                                   'diff_axis/diameter', 'diff_intensity']]

# create a list of rejected cells
props_rejected.set_index(['imageID', 'cell'], inplace=True)
rejectedCells = props_rejected.index.values.tolist()

# new trajectories df with rejected cells removed
raw_trajectories.set_index(['imageID', 'cell'], inplace=True)
clean_trajectories = raw_trajectories.drop(rejectedCells, errors='ignore')

# look if it worked
# some cells have properties but no trajectories detected so it's expected to have less then 577 cells removed from df
len(raw_trajectories.index.unique()) - len(clean_trajectories.index.unique())


## FILTER THE TRAJECTORIES # try to find short spurious trajectories seen in the images

fig, ax = plt.subplots(figsize=(10, 5))
ax = sns.scatterplot(x='time', y='distance total',
                     data=clean_trajectories.loc[clean_trajectories['frame']==clean_trajectories['frame max']],
                     alpha=0.2)

plt.show()

# try different lifetime thresholds
fig, ax = plt.subplots(ncols=3, nrows=1,
                       figsize=(20, 5))

for n, i in enumerate(['time', 'distance total', 'distance effective from centroid']):
    sns.distplot(clean_trajectories.loc[clean_trajectories['frame'] == clean_trajectories['frame max']][i],
                 ax=ax[n], hist=False, color='yellow', label='no')

    sns.distplot(clean_trajectories.loc[(clean_trajectories['frame'] == clean_trajectories['frame max'])
                                        & (clean_trajectories['time'] >= 6)][i],
                 ax=ax[n], hist=False, color='gold', label='>=6')

    sns.distplot(clean_trajectories.loc[(clean_trajectories['frame'] == clean_trajectories['frame max'])
                                        & (clean_trajectories['time'] >= 8)][i],
                 ax=ax[n], hist=False, color='orange', label='>=8')

    sns.distplot(clean_trajectories.loc[(clean_trajectories['frame'] == clean_trajectories['frame max'])
                                        & (clean_trajectories['time'] >= 10)][i],
                 ax=ax[n], hist=False, color='coral', label='>=10')

    sns.distplot(clean_trajectories.loc[(clean_trajectories['frame'] == clean_trajectories['frame max'])
                                        & (clean_trajectories['time'] >= 12)][i],
                 ax=ax[n], hist=False, color='crimson', label='>=12')

    sns.distplot(clean_trajectories.loc[(clean_trajectories['frame'] == clean_trajectories['frame max'])
                                        & (clean_trajectories['time'] >= 14)][i],
                 ax=ax[n], hist=False, color='mediumvioletred', label='>=14')

    sns.distplot(clean_trajectories.loc[(clean_trajectories['frame'] == clean_trajectories['frame max'])
                                        & (clean_trajectories['time'] >= 16)][i],
                 ax=ax[n], hist=False, color='darkviolet', label='>=16')

    sns.distplot(clean_trajectories.loc[(clean_trajectories['frame'] == clean_trajectories['frame max'])
                                        & (clean_trajectories['time'] >= 18)][i],
                 ax=ax[n], hist=False, color='mediumblue', label='>=18')

plt.show()
# thresholds under 10 seconds seem to affect the distribution of distances, it fits with the images
# filter trajectories that last less than 10s
clean_trajectories = clean_trajectories.loc[(clean_trajectories['lifetime']) >= 10]

# check if it worked
clean_trajectories['lifetime'].min()


## COUNT DETECTED TRAJECTORIES PER STRAIN

# count cells
Ncells = clean_trajectories.loc[~clean_trajectories.index.duplicated(keep='first')].groupby('strain')['particle'].count().to_frame('cells')

# count total number of particles detected
Ntotal_particles = clean_trajectories[clean_trajectories['time'] == 0] # subset first frame of each trajectories
Ntotal_particles = Ntotal_particles.groupby('strain')['particle'].count().to_frame('total particles')

# count number of particles with lifetime >=180
Nlifetime180 = clean_trajectories[clean_trajectories['time'] == 180]
Nlifetime180 = Nlifetime180.groupby('strain')['particle'].count().to_frame('lifetime >= 180')

# count number of new dots appearing after the first frame
Nnew_particles = clean_trajectories[(clean_trajectories['frame min'] != 0)
                                   & (clean_trajectories['time'] == 0)] # first frame of new trajectories
Nnew_particles = Nnew_particles.groupby('strain')['particle'].count().to_frame('new particles')

# count initial dots that disappears
Ndisappearing = clean_trajectories[(clean_trajectories['frame min'] == 0)
                                     & (clean_trajectories['frame max'] != 180)
                                     & (clean_trajectories['time'] == 0)]
Ndisappearing = Ndisappearing.groupby('strain')['particle'].count().to_frame('initial particles that disappear')

# new df with only complete trajectories
complete_trajectories = clean_trajectories[(clean_trajectories['frame min'] != 0)
                                           & (clean_trajectories['frame max'] != 180)] # remove trajectories present at the beginning or at the end

# new df with last values of each complete trajectories
lastframe = complete_trajectories.loc[complete_trajectories['frame'] == complete_trajectories['frame max']]

# count number of complete trajectories
Ncomplete = lastframe.groupby('strain')['particle'].count().to_frame('complete trajectories')

# nb of cells and trajectories for each strain
counts = pd.concat([Ncells, Ntotal_particles, Nlifetime180, Nnew_particles, Ndisappearing, Ncomplete], axis=1, sort=False)

# calculate nb per cell
counts['180/cell'] = counts['lifetime >= 180'] / counts['cells']
counts['new/cell'] = counts['new particles'] / counts['cells']
counts['start/cell'] = counts['initial particles that disappear'] / counts['cells']
counts['comp/cell'] = counts['complete trajectories'] / counts['cells']
counts['total/cell'] = counts['total particles'] / counts['cells']
counts['incomp/cell'] = (counts['total particles'] - counts['complete trajectories']) / counts['cells']

# save dataframes
clean_trajectories.to_csv('SLA1GFPclean_trajectories.csv')
complete_trajectories.to_csv('SLA1GFPcomplete_trajectories.csv')
lastframe.to_csv('SLA1GFPlastframe.csv')
counts.to_csv('SLA1GFPcounts.csv')


## STACKED PLOTS COMPLETE VS INCOMPLETE TRAJECTORIES PER STRAIN (Figures 4H and S5G)

# save csv
stacked = counts.loc[:, ['comp/cell', 'incomp/cell']]
stacked.rename(columns={'comp/cell':'complete', 'incomp/cell':'incomplete'}, inplace=True)
stacked.to_csv('SLA1GFPCompleteAndIncomplete.csv')

# stuffed strains (Figure S5G)
stuffed = ['1|2|3', 'D|2|3', '1|D|3', '1|2|D', 'D|D|D']
ind = np.arange(5)
sub = stacked.loc[stuffed, :]

fig = plt.figure(figsize=(2.8, 5))
plt.bar(ind, sub['complete'], width=0.6, bottom=sub['incomplete'], label='Complete', color='silver')
plt.bar(ind, sub['incomplete'], width=0.6, label='Incomplete', color='grey')
plt.ylabel('Trajectories per cell')
plt.xticks(ind, stuffed, rotation=90)
plt.legend(loc="upper right")

fig.savefig('plots/nTrajectories-stuffed.svg', bbox_inches='tight')

# shuffles (Figure 4H)
shuffles = ['1|2|3', 'D|D|D', '2|1|3', '1|3|2', '3|2|1', '2|3|1', '3|1|2']
ind = np.arange(7)
sub = stacked.loc[shuffles, :]

fig = plt.figure(figsize=(4, 5))
plt.bar(ind, sub['complete'], width=0.6, bottom=sub['incomplete'], label='Complete', color='silver')
plt.bar(ind, sub['incomplete'], width=0.6, label='Incomplete', color='grey')
plt.ylabel('Trajectories per cell')
plt.xticks(ind, shuffles, rotation=90)
plt.legend(loc="upper right")

fig.savefig('plots/nTrajectories-shuffles.svg', bbox_inches='tight')


## LOOK AT DIFFERENCE BETWEEN STRAINS

# calculate speed and linearity
lastframe['effective speed'] = lastframe['distance effective from centroid']/lastframe['lifetime']
lastframe['linearity'] = lastframe['distance effective']/lastframe['distance total']

# make violin plots
fig, ax = plt.subplots(ncols=3, nrows=2,
                       figsize=(20, 15),
                       sharey='row')
plt.subplots_adjust(wspace=0.1, hspace=0.2)

order = ['1|2|3',
       'D|2|3', '2|2|3', '3|2|3',
       '1|D|3', '1|3|3',
       '1|2|D', '1|2|1', '1|2|2',
       'D|D|3', '2|1|3',
       'D|2|D', '3|2|1',
       '1|D|D', '1|3|2',
       'D|D|D', '2|3|1', '3|1|2']

# for n, i in enumerate(['distance effective', 'distance total', 'distance effective from centroid',
#                        'lifetime', 'effective speed', 'linearity']):
sns.violinplot(y='strain', x='distance effective', data=lastframe, ax=ax[0,0], order=order, showfliers=False)
sns.violinplot(y='strain', x='distance total', data=lastframe, ax=ax[0,1], order=order, showfliers=False)
sns.violinplot(y='strain', x='distance effective from centroid', data=lastframe, ax=ax[0,2], order=order, showfliers=False)
sns.violinplot(y='strain', x='lifetime', data=lastframe, ax=ax[1,0], order=order, showfliers=False)
sns.violinplot(y='strain', x='effective speed', data=lastframe, ax=ax[1,1], order=order, showfliers=False)
sns.violinplot(y='strain', x='linearity', data=lastframe, ax=ax[1,2], order=order, showfliers=False)

plt.show()

# save csv for boxplots (Figures S5H and S5I)
df_boxplots = lastframe[['strain','lifetime', 'linearity']]
df_boxplots.to_csv('SLA1GFPLifetimeAndLinearity.csv')


## PLOT MOVEMENTS TOWARDS CENTROID OVER TIME (Figures 4G and S5F)

# calculate effective distance towards centroid
complete_trajectories.reset_index(inplace=True)
complete_trajectories.set_index(['imageID', 'cell', 'particle'], inplace=True)
complete_trajectories['initial distance from centroid'] = 0.0
complete_trajectories.update(complete_trajectories.loc[complete_trajectories['time'] == 0,
                                                       'distance effective from centroid per frame'
                                                      ].to_frame('initial distance from centroid'))
complete_trajectories['distance effective towards centroid'] = complete_trajectories['initial distance from centroid'] \
                                                               - complete_trajectories['distance effective from centroid per frame']

df_lineplots = complete_trajectories.reset_index().groupby(['strain', 'time'])['distance effective towards centroid']\
                                                  .agg(['mean', 'count'])
df_lineplots.reset_index(inplace=True)

# save csv
df_lineplots.rename(columns={'mean':'mean effective distance towards centroid', 'count':'n'}, inplace=True)
df_lineplots.to_csv(f'SLA1GFPDistances.csv')

# plot stuffed strains (Figure S5F)
fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(7, 5))
plt.axis([0, 120, -50, 100])
ax[0, 0].set_xlabel('Time (sec)')
ax[0, 0].set_ylabel('Effective distance towards centroid (nm)')
ax[0, 0].axhline(0, color='grey', linewidth=0.5, linestyle='-.')
colors = ['black', 'red', 'purple', 'blue', 'green']
n_dots = [] # to save initial number of dots
dots95 = pd.DataFrame() # to save coordinates where 95% of dots disappeared

for n, i in enumerate(stuffed):
    subset = df_lineplots[df_lineplots['strain'] == i]
    x = subset['time']
    y = subset['mean effective distance towards centroid']
    subset['count_proportion'] = subset['n'] / subset['n'].max()
    count = subset['count_proportion']
    subset['count_dist95'] = abs(subset['count_proportion'] - 0.05)

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # map segments
    lc = LineCollection(segments, linewidths=3,
                        cmap=LinearSegmentedColormap.from_list('myCmap', ['white', colors[n]]))

    # set values used for colormapping
    lc.set_array(count)
    line = ax[0, 0].add_collection(lc)

    # add black colorbar
    if n == 0:
        fig.colorbar(line, use_gridspec=False, anchor=(0., 0.), label='Proportion of particles remaining', shrink=0.45)

    n_dots.append(subset['n'].max()) # initial number of dots
    index95 = subset.loc[subset['count_dist95'].idxmin()]

    # add "+" where 95% of dots disappeared (coordinates are normalised for plt.axis([0, 120, -50, 100]))
    ymin = ((index95['mean effective distance towards centroid'] + 48.5) / 150)
    ymax = ((index95['mean effective distance towards centroid'] + 51.5) / 150)
    xmin = ((index95['time'] - 1) / 120)
    xmax = ((index95['time'] + 1) / 120)
    ax[0, 0].axvline(index95['time'],
                     ymin=ymin, ymax=ymax,
                     color=colors[n], linewidth=1)
    ax[0, 0].axhline(index95['mean effective distance towards centroid'],
                     xmin=xmin, xmax=xmax,
                     color=colors[n], linewidth=1)

# make legend
name_to_color = {f'{stuffed[0]}  (n = {n_dots[0]})': colors[0],
                 f'{stuffed[1]}  (n = {n_dots[1]})': colors[1],
                 f'{stuffed[2]}  (n = {n_dots[2]})': colors[2],
                 f'{stuffed[3]}  (n = {n_dots[3]})': colors[3],
                 f'{stuffed[4]}  (n = {n_dots[4]})': colors[4]}
patches = [Patch(color=v, label=k) for k, v in name_to_color.items()]
plt.legend(handles=patches, bbox_to_anchor=(1.45, 1))

fig.savefig(f'plots/Distances-stuffed.svg', bbox_inches='tight')

# plot shuffles (Figure 4G)
fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(7, 5))
plt.axis([0, 120, -50, 100])
ax[0, 0].set_xlabel('Time (sec)')
ax[0, 0].set_ylabel('Effective distance towards centroid (nm)')
ax[0, 0].axhline(0, color='grey', linewidth=0.5, linestyle='-.')
colors = ['black', 'red', 'purple', 'blue', 'darkcyan', 'green', 'darkorange']
n_dots = [] # to save initial number of dots
dots95 = pd.DataFrame() # to save coordinates where 95% of dots disappeared

for n, i in enumerate(shuffles):
    subset = df_lineplots[df_lineplots['strain'] == i]
    x = subset['time']
    y = subset['mean effective distance towards centroid']
    subset['count_proportion'] = subset['n'] / subset['n'].max()
    count = subset['count_proportion']
    subset['count_dist95'] = abs(subset['count_proportion'] - 0.05)

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # map segments
    lc = LineCollection(segments, linewidths=3,
                        cmap=LinearSegmentedColormap.from_list('myCmap', ['white', colors[n]]))

    # set values used for colormapping
    lc.set_array(count)
    line = ax[0, 0].add_collection(lc)

    # add black colorbar
    if n == 0:
        fig.colorbar(line, use_gridspec=False, anchor=(0., 0.), label='Proportion of particles remaining', shrink=0.45)

    n_dots.append(subset['n'].max()) # initial number of dots
    index95 = subset.loc[subset['count_dist95'].idxmin()]

    # add "+" where 95% of dots disappeared (coordinates are normalised for plt.axis([0, 120, -50, 100]))
    ymin = ((index95['mean effective distance towards centroid'] + 48.5) / 150)
    ymax = ((index95['mean effective distance towards centroid'] + 51.5) / 150)
    xmin = ((index95['time'] - 1) / 120)
    xmax = ((index95['time'] + 1) / 120)
    ax[0, 0].axvline(index95['time'],
                     ymin=ymin, ymax=ymax,
                     color=colors[n], linewidth=1)
    ax[0, 0].axhline(index95['mean effective distance towards centroid'],
                     xmin=xmin, xmax=xmax,
                     color=colors[n], linewidth=1)

# make legend
name_to_color = {f'{shuffles[0]}  (n = {n_dots[0]})': colors[0],
                 f'{shuffles[1]}  (n = {n_dots[1]})': colors[1],
                 f'{shuffles[2]}  (n = {n_dots[2]})': colors[2],
                 f'{shuffles[3]}  (n = {n_dots[3]})': colors[3],
                 f'{shuffles[4]}  (n = {n_dots[4]})': colors[4],
                 f'{shuffles[5]}  (n = {n_dots[5]})': colors[5],
                 f'{shuffles[6]}  (n = {n_dots[6]})': colors[6]}
patches = [Patch(color=v, label=k) for k, v in name_to_color.items()]
plt.legend(handles=patches, bbox_to_anchor=(1.45, 1))

fig.savefig(f'plots/Distances-shuffles.svg', bbox_inches='tight')


## PLOT EXAMPLES OF TRAJECTORIES IN CELLS (Figure S5E)

# make function to plot the trajectories of the particles for a single cell
def plot_trajectories(ID, cell):
    sub = complete_trajectories.loc[(complete_trajectories['imageID'] == ID)
                                    & (complete_trajectories['cell'] == cell)]

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.set_title(sub['strain'].min())
    ax.set_xlim(0, 150)
    ax.set_ylim(0, 150)
    ax.set(xticklabels=[])
    ax.set(yticklabels=[])

    # load cell perimeter
    img_bright = np.load(f'AnalysedImages/{ID}/cells/cell{cell:08d}//cellbright.npy')
    ax.contour(img_bright, [0.5], linewidths=1, linestyles='dashed', colors='grey')

    # add trajectories
    for i in sub.particle.unique():
        x = sub[sub['particle'] == i]['x']
        y = sub[sub['particle'] == i]['y']
        time = sub[sub['particle'] == i]['time']
        norm = plt.Normalize(sub[sub['particle'] == i]['time'].min(),
                             sub[sub['particle'] == i]['time'].max())

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, linewidths=2, cmap='viridis', norm=norm)
        lc.set_array(time)
        line = ax.add_collection(lc)
        plt.axis([0, 150, 0, 150])

    # add colorbar
    fig.colorbar(line, shrink=0.6, label='Time (sec)')

    fig.savefig(f'plots/Trajectories-{ID}_c{cell}.svg', bbox_inches='tight')

# plot trajectories over time (Figure S5E)
complete_trajectories.reset_index(inplace=True)
plot_trajectories('WTWTWT_1211_006', 10)
