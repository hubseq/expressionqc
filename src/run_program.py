import os, subprocess, sys
import pandas as pd
import scipy as sp
import numpy as np
from sklearn.decomposition import PCA
import plotly.express as px
sys.path.append('/global_utils/src/')
# sys.path.append('global_utils/src/')
import module_utils
from sklearn.preprocessing import StandardScaler
from scipy import stats
import plotly.figure_factory as ff
import matplotlib.pyplot as plt
from io import BytesIO
import base64

def array2matrix( array1d, dim ):
    """ Converts 1D array of values to a 2D matrix of size (dim, dim). 
    For plotting correlation heatmap.
    """
    array2d = []
    dim_og = dim
    while dim > 0:
        array2d.append([0]*(dim_og-dim) + array1d[0:dim])
        dim = dim - 1
        array1d = array1d[dim:]
    return array2d

def remove_low_counts( x, y, t = 10 ):
    """ Remove all counts below a threshold (default 10)
    """
    x_new, y_new = [], []
    for i in range(len(x)):
        if not((x[i] < t and y[i] < t)): # or (x[i] < t and y[i] < 10*t) or (x[i] < 10*t and y[i] < t)):
            x_new.append(x[i])
            y_new.append(y[i])
    return x_new, y_new

def correlation_heatmap( x, s, cscale, samples ):
    d = dict(zip(list(range(0,len(samples))), samples))
    
    R_matrix = []
    for i in range(len(s)):
        R_row = []
        for j in range(len(s)):
            if j >= i:
                x_new, y_new = remove_low_counts(list(x.iloc[s[i],:]),list(x.iloc[s[j],:]))
                R, P = stats.spearmanr(x_new, y_new)
                R_row.append(R)
            else:
                R_row.append(0)
        R_matrix.append(R_row)
    
    labels = list(map(lambda x: d[x], s))
    
    fig = ff.create_annotated_heatmap(
        np.array(R_matrix),
        x = samples,
        y = samples,
        annotation_text = np.around(R_matrix, decimals=2),
        hoverinfo='z',
        colorscale=cscale
        )
    fig.update_layout(
        autosize=False,
        width=600,
        height=600)
    fig.write_image('expressionqc.correlation_heatmap.png')
    
    return R_matrix, fig


def scatter_plots( x, samples, s, group_name ):
    d = dict(zip(list(range(0,len(samples))), samples))
    
    R_all = []
    print('SCATTER_PLOTSS: {} {} {} {}'.format(str(x), str(samples), str(s), str(group_name)))
    fig, axs = plt.subplots(len(s), len(s), figsize=(20, 20))
    plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.4)
    plt.title(group_name)
    for i in range(len(s)):
        for j in range(len(s)):
            if j >= i:
                x_new, y_new = remove_low_counts(list(x.iloc[s[i],:]),list(x.iloc[s[j],:]))
                R, P = stats.spearmanr(x_new, y_new)
                axs[i][j].set_xscale('log')
                axs[i][j].set_yscale('log')
                axs[i][j].set_xlabel(d[s[i]])
                axs[i][j].set_ylabel(d[s[j]])
                axs[i][j].scatter(list(x.iloc[s[i],:]),list(x.iloc[s[j],:]))
                axs[i][j].text(1, 10**5, 'R={}'.format(str(round(R,3))))
                if j != i:
                    R_all.append(R)
    plt.savefig('expressionqc.scatter.{}.png'.format(group_name))
    return np.mean(R_all), plt


def scatter_plots_pairwise( x, samples, s1, s2, group_name1, group_name2 ):
    d = dict(zip(list(range(0,len(samples))), samples))    
    
    R_all = []
    if len(s1) > 1 or len(s2) > 1:
        fig, axs = plt.subplots(len(s1), len(s2), figsize=(20, 20))
        plt.subplots_adjust(left=0.1,
                            bottom=0.1, 
                            right=0.9, 
                            top=0.9,
                            wspace=0.4,
                            hspace=0.4)
        plt.title('{} vs {}'.format(group_name1, group_name2))
        for i in range(len(s1)):
            for j in range(len(s2)):
                x_new, y_new = remove_low_counts(list(x.iloc[s1[i],:]),list(x.iloc[s2[j],:]))
                R, P = stats.spearmanr(x_new, y_new)
                axs[i][j].set_xscale('log')
                axs[i][j].set_yscale('log')
                axs[i][j].set_xlabel(d[s1[i]])
                axs[i][j].set_ylabel(d[s2[j]])
                axs[i][j].scatter(list(x.iloc[s1[i],:]),list(x.iloc[s2[j],:]))
                axs[i][j].text(1, 10**5, 'R={}'.format(str(round(R,3))))
                if j != i:
                    R_all.append(R)
    else:
        # only 1 comparison
        plt.title('{} vs {}'.format(group_name1, group_name2))
        x_new, y_new = remove_low_counts(list(x.iloc[s1[0],:]),list(x.iloc[s2[0],:]))
        R, P = stats.spearmanr(x_new, y_new)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(d[s1[0]])
        plt.ylabel(d[s2[0]])
        plt.scatter(list(x.iloc[s1[0],:]),list(x.iloc[s2[0],:]))
        plt.text(1, 10**5, 'R={}'.format(str(round(R,3))))
        R_all.append(R)
    plt.savefig('expressionqc.scatter.{}-vs-{}.png'.format(group_name1, group_name2))
    return np.mean(R_all), plt


def ReadsPerGene2Matrix( count_files, count_file_type='STAR', experiment_type='unstranded'):
    """
    count_file_type: 'STAR', 'salmon', 'featurecounts'
    	which aligner or program was used to generate counts files.

    return_type:
    column = each sample has counts in a separate column, like this:
    	gene,sample1,sample2,sample3,...
    	ENSG001,5,8,12,6
    	ENSG002,20,15,3,8
    	...
    tallskinny = stacked, like this:
    	sample,gene,count
    	sample1,ENSG001,5
    	sample1,ENSG002,20
    	sample2,ENSG001,8
    	sample2,ENSG002,15
    	...
    transpose = samples are the rows and genes are the columns, like this:
    	sample,ENSG001,ENSG002,...
    	sample1,5,20,...
    	sample2,8,15,...
    
    NOTE: all types are output as files, but only one is returned for PCA analysis.
    """
    if count_file_type == 'STAR':
        if experiment_type=='unstranded':
            WHICH_COL = 1
        elif experiment_type=='stranded':
            WHICH_COL = 2
        elif experiment_type=='reverse-stranded':
            WHICH_COL = 3
        
    for output_type in ['column','transpose','stacked']:
        counts_matrix_dict = {}
        genes = []
        samples = []
        firstFile = True        
        for cf in count_files:
            print('Parsing {}'.format(str(cf)))
            if count_file_type == 'STAR':
                with open(cf,'r') as f:
                    samplename = str(cf.split('/')[-1]).split('.')[0]
                    samples.append(samplename)
                    for r in f:
                        rt = r.rstrip(' \t\n').split('\t')
                        if r[0:2] != 'N_':
                            if output_type == 'column':
                                if firstFile == True:
                                    genes.append(rt[0])
                                    counts_matrix_dict[rt[0]] = [rt[WHICH_COL]]  # gene: counts for each sample
                                else:
                                    counts_matrix_dict[rt[0]] += [rt[WHICH_COL]]
                            elif output_type == 'transpose' or output_type in ['tallskinny', 'stacked']:
                                if firstFile == True:
                                    genes.append(rt[0])
                                # sample: counts for each gene                                    
                                if samplename not in counts_matrix_dict:
                                    counts_matrix_dict[samplename] = [rt[WHICH_COL]]  
                                else:
                                    counts_matrix_dict[samplename] += [rt[WHICH_COL]]
                firstFile = False
            
        with open('expressionqc.counts_matrix.{}.csv'.format(output_type),'w') as fout:
            if output_type == 'column':
                fout.write('gene,{}\n'.format(','.join(samples)))
                for gene, counts in counts_matrix_dict.items():
                    fout.write('{},{}\n'.format(gene, ','.join(counts)))
            elif output_type == 'transpose':
                fout.write('sample,{}\n'.format(','.join(genes)))
                for sample, counts in counts_matrix_dict.items():
                    fout.write('{},{}\n'.format(sample,','.join(counts)))
            elif output_type in ['tallskinny', 'stacked']:
                fout.write('sample,gene,count\n')
                for sample, counts in counts_matrix_dict.items():
                    for i in range(0,len(counts)):
                        fout.write('{},{},{}\n'.format(sample,genes[i],str(counts[i])))
    
    return 'expressionqc.counts_matrix.transpose.csv', samples


def runPCA( counts_json, samples, groups ):
    def groups2colors( g ):
        """ Given list of group names, assigns colors
        """
        c = ['#2E91E5','#E15F99','#1CA71C','#222A2A','#FBE426', '#16FF32', \
             '#B68100', '#750D86', '#EB663B', '#B2828D', '#511CFB', '#00A08B', '#AF0038']
        which_c = -1
        g_previous = ''
        c_out = []
        for g_current in g:
            if g_current != g_previous:
                which_c = which_c + 1 if which_c < len(c) - 1 else 0
            c_out.append(c[which_c])
            g_previous = g_current
        return c_out

    counts = counts_json['counts']
    target = counts_json['target']
    features = counts_json['features']
    print(counts)
    
    # perform standard scalar normalization to normalize our feature set
    # Standardize the features
    counts_norm = StandardScaler().fit_transform(counts)
    
    ########## Instantiate PCA 2D ##########
    pca = PCA(n_components=2)
    
    # Fit PCA to features
    principalComponents = pca.fit_transform(counts_norm)
    
    ##### explained variance ratio #####
    explained_variance = pca.explained_variance_ratio_
    with open('expression.counts_pca_explained_variance.csv','w') as fout:
        fout.write('PC1_perc_variance,PC2_perc_variance\n')
        fout.write(','.join(list(map(str,explained_variance)))+'\n')
    
    ##### Principal Components DataFrame #####
    # Create a new dataset from principal components 
    df_pc2d = pd.DataFrame(data = principalComponents, 
                      columns = ['PC1', 'PC2'])  # 'PC3'
    
    principal_df = pd.concat([df_pc2d, target], axis=1)
    principal_df['Group'] = groups
    principal_df.to_csv('expressionqc.counts_pca.csv')
    
    # Output data points to CSV
    counts.to_csv('expressionqc.counts.csv')
    
    # Create PCA plot
    colors = groups2colors(groups)
    fig = px.scatter(principal_df, x='PC1', y='PC2', color="Group", text=samples, labels={"Group": "Group"})
    fig.update_traces(textposition='top center', marker_size=10)
    fig.write_image("counts_pca.png")

    ########## Instantiate PCA 3D ########## 
    pca_3d = PCA(n_components=3)
    
    # Fit PCA to features
    principalComponents_3d = pca_3d.fit_transform(counts_norm)
    # Create a new dataset from principal components 
    df_3d = pd.DataFrame(data = principalComponents_3d, 
                         columns = ['PC1', 'PC2', 'PC3'])

    principal_df_3d = pd.concat([df_3d, target], axis=1)
    principal_df_3d['Group'] = groups
    
    # Explained variance %
    # explained variance ratio
    explained_variance_3d = pca_3d.explained_variance_ratio_
    with open('expression.counts_pca_3d_explained_variance.csv','w') as fout:
        fout.write('PC1_perc_variance,PC2_perc_variance,PC3_perc_variance\n')
        fout.write(','.join(list(map(str,explained_variance_3d)))+'\n')
        
    # Visualization
    fig2 = px.scatter_3d(principal_df_3d, x='PC1', y='PC2', z='PC3', color="Group", text=samples, labels={"Group": "Group"})
    fig2.update_layout(
        autosize=False,
        width=800,
        height=1200)
    fig2.update_traces(textposition='top center', marker_size=10)
    fig2.write_image("counts_pca_3d.png")
    
    return fig, fig2


def createScatterPlots( counts, samples, groups ):
    # scatter plot figure objects
    plotfigs = []
    # get within group comparisons
    within_group_comps = []
    new_group_list = []
    prev_group = ''
    for i in range(0,len(groups)):
        if (groups[i] != prev_group and new_group_list != []):
            within_group_comps.append(new_group_list)
            new_group_list = [i]
            prev_group = groups[i]
        else:
            new_group_list.append(i)
            prev_group = groups[i]
        
        # append the last group
        if i >= len(groups)-1:
            within_group_comps.append(new_group_list)

    print('WITHIN GROUP COMPS: {}'.format(str(within_group_comps)))
    # within group scatter plots
    for within_group_indexes in within_group_comps:
        if len(within_group_indexes) > 1:
            R_avg, plotfig = scatter_plots( counts, samples, within_group_indexes, groups[within_group_indexes[0]] )
            plotfigs.append([plotfig, 'matplotplib'])
    
    # cross group scatter plots
    for i in range(0,len(within_group_comps)):
        for j in range(i,len(within_group_comps)):
            if j > i:
                R_avg, plotfig = scatter_plots_pairwise( counts, samples, within_group_comps[i], within_group_comps[j], \
                                                groups[within_group_comps[i][0]], groups[within_group_comps[j][0]] )
                plotfigs.append([plotfig, 'matplotlib'])
    
    # correlation heatmap
    R_matrix, plotfig = correlation_heatmap( counts, list(range(0,len(samples))), 'hot', samples )
    plotfigs.append([plotfig, 'plotly'])

    return plotfigs


def outputSampleGroupMetadata( samples, groups ):
    """ Creates a CSV containing samples and their groups
    """
    if type(samples)!=type([]) or type(groups)!=type([]) or len(samples) != len(groups):
        out_name = ''
        print('ERROR in outputSampleGroupMetadata(): samples and groups lists not the same length.')
        print('SAMPLES: {}'.format(str(samples)))
        print('GROUPS: {}'.format(str(groups)))
    else:
        out_name = 'expressionqc.samplegroups.csv'
        with open(out_name,'w') as fout:
            fout.write('sample_id,group\n')
            for i in range(0,len(samples)):
                fout.write('{},{}\n'.format(str(samples[i]), str(groups[i])))
    return out_name

def add_image_to_html( p, p_type ):
    """ Given <plot_figure> and <plot_type>, returns an img src tag with embedded image
    """
    imgfile = BytesIO()
    if p_type == 'plotly':
        p.write_image(imgfile, format='png') # matplotlib is savefig
    elif p_type == 'matplotlib':
        p.savefig(imgfile)        
    encoded = base64.b64encode(imgfile.getvalue()).decode('utf-8')
    img_tag = '<img src=\'data:image/png;base64,{}\'>'.format(encoded) + '<br>'
    return img_tag


def plots_to_html( plots_object_list, html_outname ):
    """ Given a list of plot objects (plotly or matplotlib), outputs these plots to an HTML
    """
    img_tags = ''
    with open(html_outname,'w') as f:
        for p_tuple in plots_object_list:
            img_tags += add_image_to_html( p_tuple[0], p_tuple[1] )
    html = '<html><body>' + img_tags + '</body></html>'
    with open( html_outname, 'w' ) as f:
        f.write(html)
    return


def run_expressionqc( arg_list ):
    """
    Parameters:
    -i <input_count_files>
    -groups <groups>
    -o <output_dir>
    -type <STAR or featurecounts or salmon>
    -name <base name for plots> : default expressionqc
    """
    print('ARG LIST: {}'.format(str(arg_list)))
    input_args = module_utils.getArgument( arg_list, '-i', 'list' )
    output_dir = module_utils.getArgument( arg_list, '-o', 'list' )
    groups = module_utils.getArgument( arg_list, '-groups', 'implicit' )
    count_file_type = module_utils.getArgument( arg_list, '-type' )
    base_name = module_utils.getArgument( arg_list, '-name', 'implicit', 'expressionqc' )
    print('here1')
    if output_dir != []:
        os.chdir( output_dir[0] )
    print('here2')        
    if count_file_type == '':
        count_file_type = 'STAR' # default
    print('here3')        
    if groups == [] or groups == '':
        # default: group1, group2, group3, etc...
        groups = []
        for i in range(1,len(input_args)+1):
            groups.append('group{}'.format(str(i)))
    elif type(groups) == str and ',' in groups:
        groups = groups.replace(' ','').split(',')
        
    print('here4: GROUPS {}'.format(str(groups)))
    if input_args != []:
        matrix_file, samples = ReadsPerGene2Matrix( input_args, count_file_type, 'unstranded' )

        # output samples and groups as a metadata file
        outputSampleGroupMetadata( samples, groups )
        
        ### Create counts matrix ###
        print('here5: SAMPLES {}'.format(str(samples)))
        df_counts_matrix = pd.read_csv(matrix_file)
        print(df_counts_matrix)
    
        # remove any read counts less than 5
        df_counts_matrix = df_counts_matrix.loc[:, (~df_counts_matrix.isin([0,1,2,3,4])).any(axis=0)]
        
        # get target, features and counts for PCA
        target = df_counts_matrix.iloc[:,0]
        counts = df_counts_matrix.iloc[:,1:]
        features = list(df_counts_matrix.columns)[1:]

        # run PCA
        plot_pca2d, plot_pca3d = runPCA( {'counts': counts, 'target': target, 'features': features}, samples, groups )
        
        print('here6')

        # create scatter plots
        plot_scatterplots_list = createScatterPlots( counts, samples, groups )

        # output all plots to HTML
        plots_to_html( [[plot_pca2d, 'plotly']] + [[plot_pca3d, 'plotly']] + plot_scatterplots_list, \
                       os.path.join(output_dir[0], base_name + '.plots.html'))
    return

if __name__ == '__main__':
    print(' even here?')
    run_expressionqc( sys.argv[1:] )
