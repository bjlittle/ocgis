import numpy as np
import itertools
from ocgis.calc.wrap.base import OcgFunctionTree, OcgCvArgFunction
from ocgis.calc.wrap import library
from ocgis.calc.wrap.library import SampleSize


class OcgCalculationEngine(object):
    '''
    wrapf : object : this is the wrapped function object
    grouping_idx : int[] : quickly split the data into groups
    ugrouping_idx : int[] : unique representation of grouping_idx for iteration
    weights : float[] : areal weights (not normalized)
    values : float[] : value vector corresponding to dim of grouping_idx. if
        it is a multivariate function, this should be a dictionary mapping
        variables to values that the base function will know how to interpret
    '''
    
    def __init__(self,grouping,timevec,funcs,raw=False,agg=False,
                 time_range=None,mode='raw'):
        self.raw = raw
        self.agg = agg
        self.time_range = time_range
        ## subset timevec if time_range is passed
        if self.time_range is not None:
            self.timevec = timevec[(timevec>=time_range[0])*
                                   (timevec<=time_range[1])]
        else:
            self.timevec = timevec
        ## convert solitary grouping arguments to list
        if type(grouping) == str: grouping = [grouping]
        self.grouping = grouping or ['day']
        ## always calculate the sample size. do a copy so functions list cannot
        ## grow in memory. only a problem when testing.
#        funcs_copy = copy(funcs)
#        funcs_copy.insert(0,{'func':'n'})
        self.funcs = self.set_funcs(funcs)
        self.funcs = funcs
        ## get the time groups
        self.dgroups,self.dtime = self.get_distinct_groups()
        ## select which value data to pull based on raw and agg arguments
        if self.raw:
            self.use_agg = False
        elif self.raw is False and self.agg is True:
            self.use_agg = True
        else:
            self.use_agg = False
        
    def set_funcs(self,funcs):
        potentials = OcgFunctionTree.get_potentials()
        for f in funcs:
            for p in potentials:
                if p[0] == f['func']:
                    f['ref'] = getattr(library,p[1])
                    break
            if 'name' not in f:
                f['name'] = f['func']
            if 'kwds' not in f:
                f['kwds'] = {}
        return(funcs)
        
    def get_distinct_groups(self):
        ## holds date components
        dparts = {'year':[],'month':[],'day':[],'idx':[]}
        ## pull date parts from date objects and append to date part dictionary
        for ii,dt in enumerate(self.timevec):
            for grp in self.grouping:
                dparts[grp].append(getattr(dt,grp))
            dparts['idx'].append(ii)
        ## convert to numpy arrays
        for key in dparts.keys(): dparts[key] = np.array(dparts[key],dtype=int)
        ## replace empty list with a list containing NoneType for nested
        ## iterator and find unique combinations.
        duni = {}
        for key in dparts.keys():
            if key is 'idx':
                continue
            elif len(dparts[key]) is 0:
                duni[key] = np.array([None])
            else:
                duni[key] = np.unique(dparts[key]).astype(int)
                
        ## make the unique group combinations
        
        ## will hold idx to groups
        dgroups = []
        dtime = {'tgid':[],'month':[],'day':[],'year':[]}
        ## the default select all array
        bidx = np.ones(len(dparts['idx']),dtype=bool)
        ## loop for combinations
        tgid = 1
        for year,month,day in itertools.product(duni['year'],duni['month'],duni['day']):
            ## idx arrays that when combined provide a group set
            check = dict(zip(['year','month','day'],[year,month,day]))
            yidx,midx,didx = [self._get_date_idx_(bidx,dparts,part,value) 
                              for part,value in check.iteritems()]
            idx = yidx*midx*didx
            ## if dates are drilling down to day, it is possible to return date
            ## combinations that are unreasonable.
            if idx.sum() == 0:
                continue
            for key,value in check.iteritems():
                dtime[key].append(value)
            dgroups.append(idx)
            dtime['tgid'].append(tgid)
            tgid += 1
        return(dgroups,dtime)
            
    def _get_date_idx_(self,bidx,dparts,part,value):
        if part in self.grouping:
            idx = dparts[part] == value
        else:
            idx = bidx
        return(idx)
    
    def execute(self,coll):
        ## tell collection which data to return
        coll._use_agg = self.use_agg
        ## flag used for sample size calculation for multivariate calculations
        has_multi = False
        ## iterate over functions
        for cid,f in enumerate(self.funcs,start=1):
            ## change behavior for multivariate functions
            if issubclass(f['ref'],OcgCvArgFunction):
                has_multi = True
                ## cv-controlled multivariate functions require collecting
                ## data arrays before passing to function.
                kwds = f['kwds'].copy()
                for key in f['ref'].keys:
                    ## the name of the variable passed in the request
                    ## that should be mapped to the named argument
                    backref = kwds[key]
                    ## pull associated data
                    dref = coll._get_value_(backref)
                    ## update dict with properly reference data
                    kwds.update({key:dref})
                ## function object instance
                ref = f['ref'](agg=self.agg,groups=self.dgroups,kwds=kwds,weights=coll.weights)
                ## store calculation value
                coll.calc_multi[f['name']] = ref.calculate()
            else:
                ## perform calculation on each variable
                for var_name,values in coll._iter_items_():
                    ## instance of function object
                    ref = f['ref'](values=values,agg=self.agg,groups=self.dgroups,kwds=f['kwds'],weights=coll.weights)
                    ## calculate the values
                    calc = ref.calculate()
                    ## store the values
                    coll.variables[var_name].calc_value.update({f['name']:calc})
                    ## update calculation identifier
                    coll.variables[var_name].cid = np.append(coll.variables[var_name].cid,cid)
        ## calculate sample size for multivariate calculation
        if has_multi:
            for ii,value in enumerate(coll.variables.itervalues()):
                if ii == 0:
                    n = value.calc_value['n'].copy()
                else:
                    n += value.calc_value['n']
            coll.calc_multi['n'] = n
        
        return(coll)
