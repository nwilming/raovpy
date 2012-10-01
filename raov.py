import os
import tempfile
import copy
import numpy as np
from rpy2 import robjects
from rpy2.robjects.vectors import DataFrame

import measures as ms

'''
The raovpy module allows to compute repeated measures ANOVA with 
the help of R's aov and Anova functions. It does not  carry
out the anova calculation itself but delegates it to R. For
this reason Rpy2 must be installed, furthermore the Rpy2 
instance needs the 'car' package. If it is not installed:
    
    >>> import rpy2.interactive as r
    >>> r.importr("utils")
    >>> r.packages.utils.install_packages('car')

(only needs to be carried out once)

Please note that you can only carry out complete within-subject
anovas with this package.
'''
       

def aov(matrix, factor_names, measure_name, robj, interactions = '+'):
    '''
    Computes a repeated measures anova in R via the 'aov' command.
    
    This function uses R's aov function. It does not compute 
    Greehnhouse-Geisser and Huynh-Feldt corrections. Use lm_anova
    for this.
    
    Input:
        matrix : ndarray
            Each dimension of the matrix corresponds to one factor.
            The first dimension must be (!) the number of subjects.
            The values in the matrix are taken as the dependent variable. 
        factor_names : list
            List with names of each factor. The ordering must correspond
            to the dimensions given by matrix.shape.
        measure_name : str
            Name of the dependnent variable.
        robj : rpy2.robjects instance
        interactions : str
            '+' for no interactions
            '*' for all interactions
    '''

    robj.r('rm(list = ls(all = TRUE)) ')
    df = make_data_frame(matrix, 
            factor_names, measure = measure_name)
    robj.globalenv['df'] = DataFrame(df)
    robj.r('attach(df)')
    formula = ''
    error = ''
    for factor in factor_names:
        robj.r('%s<-factor(df$%s)'%(factor,factor))
        formula = formula + interactions + factor
        error = error + '*' + factor
    formula,error  = formula[1:], error[1:]
    formula = 'aov.out <- aov(%s  ~ %s + Error(subject/(%s), data=df))'%(measure_name, formula, error)
    robj.r(formula)
    print(robj.r('summary(aov.out)'))
    robj.r('detach(df)')


def lm_anova(fms, factor_dict, interactions = '+'):
    '''
    Calculates a repeated measure anova using R and the Anova function.

    This function constructs a linear model of the data and then uses 
    the Anova function from the 'car' package to compute the appropriate 
    statistics.
    
    Input:
        fms : dictionary
            Contains a datamat per subject.
        factor_dict : dictionary
            Contains factors as keys, and factor levels as values.
        interactions : string
            '+' for no interactions
            '*' for all interactions
    '''

    cell_names = set() 
    cell_fms = []
    for name, value in fms.iteritems():
        cells, factors, factor_names = filter_by_dict(value, factor_dict)
        if len(cell_names) == 0:
            cell_names = set(cells.keys())
        else:
            assert len(cell_names - set(cells.keys())) == 0
        f = lambda v: ms.percent_correct(ms.get_contingency_table(v))
        pc = dv_from_cells(cells, f)
        cell_fms.append(pc)
    cell_names = list(cell_names)
    datafile = tempfile.NamedTemporaryFile('w', delete = False)
    factorfile = tempfile.NamedTemporaryFile('w', delete = False)
    cells2file(cell_fms, cell_names, factor_names, factors, datafile, factorfile=factorfile)    
    datafile.close()
    factorfile.close()
    _lm_anova(datafile.name, factorfile.name, factor_names, robjects, interactions)
    os.remove(datafile.name)
    os.remove(factorfile.name)


def _lm_anova(dfile, facfile, factors, robj, interactions = '+'):
    '''
    Computes a repeated measures anova in R.

    Input:
        dfile : string
            Path to file that contains data for the anova in R.
            Data in dfile must have 'wide' format, that is one
            subject per line. Every cell of the anova is one column,
            data is cells are comma separated. The first row contains
            cell names.
        facfile : string
            Filename of the file that contains factor definitions.
            The file contains one column per factor, each row indicates
            the factor level of a cell. The cells are in the same order
            as in dfile. That is, moving along the columns in dfile 
            corresponds to moving along the rows in facfile.
        factors : list
            List of factor names.
        intteractions : string
            '+' if no interactins should be estimated, '*' if interactions
            should be estimated.
    '''
    robj.r('rm(list = ls(all = TRUE)) ')
    robj.r('design_mat <- read.csv("%s", header=T)'%facfile)
    robj.r('attach(design_mat)')
    robj.r('library(car)')
    robj.r("df <- read.csv('%s', header=T)"%dfile)
    robj.r("lmmod <- lm(as.matrix(df) ~ 1)")
    formula = factors[0]
    for f in factors[1:]:
        formula = formula + interactions + f
    print formula
    robj.r("a <- Anova(lmmod, idata=design_mat, idesign =~%s, type='III')"%formula)
    print robj.r('summary(a, multivariate=FALSE)')
    robj.r('detach(design_mat)')

def filter_by_dict(fm, d):
    '''
    Filter a datamat by combinations of fields, a.k.a cells from datamat.

    d is dict that contains as keys field names and as 
    values a list that contains admissible field values.

    The function creates the set of all possible 
    field value pairs in d and filters the fixmat once
    for every combination.

    Example:
    d = {'session':[0,3], 'name':[1,3]}
    t = filter_by_dict(fm,d)
    t = {'session0name1':datamat, # contains only values wheren session==0 and name==1
         'session0name3':datamat, # session==0 and name==3
         'session3name1':datamat, # session==3 and name==1
         'session3name3':datamat} # session==3 and name==3
    '''
    results = {}
    factors = {}
    factor_names = []
    def drill_down(fm,d, name, factor_list):
        if d == {}:
            factors[name] = factor_list
            results[name] = fm
        else:
            # Pick the first key in d
            key = np.sort(d.keys())[0]
            if not key in factor_names:
                factor_names.append(key)
            values = d[key]
            for v in values:
                m = copy.copy(d)
                del m[key]
                v2 = copy.copy(factor_list)
                v2.append(v)
                drill_down(fm[fm.field(key)==v],m, name+'%s%i'%(key,v), v2)
    
    drill_down(fm, d, '', [])
    return results, factors, factor_names

def dv_from_cells(cell_dict, function):
    '''
    Computes the dependent variable for every cell.

    Input:
        cell_dict : dictionary
            Contains as keys cell names, and as values datamats that contain 
            all data for the given cell.
    Returns:
        cell_dict : dictionary
            Contains as keys cell names, and as values the dependent measu
    '''
    cell_dict = dict((key, function(value)) for (key,value) in cell_dict.iteritems())
    return cell_dict




def cells2file(cell_fms, cell_names, factor_names, factors, datafile, factorfile):
    """
    Exports dv values to files that can be read by R and SPSS.

    Inpute:
        cell_fms : List 
            Each entry is a dictionary with cell names as keys and 
            the dependent variable as value.
        cell_names : list
            Names of the cells, e.g. the keys of the subject dictionaries in
            cell_fms.
        factor_names : list 
            Names of the factors
        factors : dictionary
            Contains for each cell a list of factor levels that 
            encodes the level of each factor for this cell. The list
            has to be ordered according to factor_names.
        filename : File object
            File to output data to. Must be writable.
        factorfile : File object
            File to output factor definitions to. Must be writable.
    """
    cell_names = np.sort(cell_names)
    # Write header with variable names
    for cell in cell_names[:-1]:
        datafile.write('%s,'%cell)
    datafile.write('%s\n'%cell_names[-1])

    # Write subject data
    for datum in cell_fms:
        for cell in cell_names[:-1]:
            datafile.write('%f,'%datum[cell])
        datafile.write('%f\n'%datum[cell_names[-1]])
    # Last thing: Create a factor file that encodes how the
    # factors change over the cells
    for f in factor_names[:-1]:
        factorfile.write('%s,'%f)
    factorfile.write('%s\n'%factor_names[-1])
    
    for cell in cell_names:
        for fv in factors[cell][:-1]:
            factorfile.write('%sv,'%fv)
        factorfile.write('%sv\n'%factors[cell][-1])


def make_data_frame(matrix, fields, measure='perc_corr'):
    '''
    Returns a dictionary that can be passed into R for data analysis.

    Data Frames will loke like this:
    Subject Field[1] Field[2] ... Field[n]
    1       x        y        ... z
    1 ...
    2       a        b        ... d
    ...

    Input:
        matrix: ndarray
            Each dimension is treated as one factor and dummy
            coded (i.e. 0:n)
        fields: iterable
            String names for each factor
        measure: string
            String Name for the dependent variable
    Output:
        Dictionary
    '''
    slice_list = [slice(d) for d in matrix.shape]
    indices = np.mgrid[slice_list]
    
    data_frame = {measure:robjects.FloatVector(matrix.flatten().tolist())}
    for index, name in zip(indices, fields):
        data_frame[name] = robjects.IntVector((1+index).flatten().tolist())

    return data_frame

