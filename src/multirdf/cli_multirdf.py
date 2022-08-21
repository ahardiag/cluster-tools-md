#!/usr/bin/env python 

"""
Run several times the save_RDF script to compute Radial Distribution Functions 

Outputs:
- png file in directory ./Figures
- csv data in directory ./Output
- log files in directory ./outputs/logs
"""

from configparser import ConfigParser
import os,re,sys
import src.multirdf.onerdf as onerdf
import argparse
from argparse import RawTextHelpFormatter
import datetime

class RunRdf(object):
    """Interface for generating parameters and running RDF script.

    Parameters
    ----------
    input_file : str, optional
        Path to input file.
        TODO: Doc intput file.
    default_parameters : dict, optional
        Default values of the job parameters.
        Example: `{"SUB": "10", "NBINS": 100}`
    """
    def __init__(self,
                 args_script,
                 default_parameters=None,
                 **kwargs
                 ):
        self.default_parameters = default_parameters
        self.args_script=args_script
        
    def read_input(self,printout=False):
        '''
        Read lines of the handle file and return a list of dictionaries
        in which each element corresponds to one set of parameters
        Removes also comments, blank lines.
        '''
        f=open(self.args_script.input_file)
        self.params_input=f.readline().rstrip('\n')
        self.params_input=re.findall(r'\w+',self.params_input)
        
        self.list_par=[]
        pattern=r'(?:\".*?\"|\S)+'
        comment_pattern=re.compile(r"^\s*#")
        blank_pattern = re.compile(r"^\s*$")

        for line in f.readlines():
            line=line.rstrip('\n')
            if comment_pattern.match(line):
                continue
            elif blank_pattern.match(line):
                continue
            else:
                params_items=re.findall(pattern,line)
                self.list_par.append({key:val \
                for (key,val) in zip(self.params_input,params_items)})
                
            if printout:
                print("Input parameters that replace default values :{}".format(self.params_input))
                print("List of dictionaries with each set of input parameters \n {}".format(self.list_par))

            assert len(self.params_input)==len(params_items), \
              'In the input file, the number of columns must match between headers and other rows !'
            
        return self.list_par
    
    def build_param(self,**kwargs):
        '''
        Build a list of objects ParamSet with all parameters to run the final script
        '''
        self.list_ps=[]
        for dic in self.list_par:
            ps=ParamSet(self.args_script)
            ps.add_par(kwargs['default_parameters'])
            ps.add_par(dic)
            self.list_ps.append(ps)
    
    def check_complete(self):
        """
        Check the necessary requirements for each ParamSet. 
        """
        for ps in self.list_ps:
            ps.check_complete()
    
    def run_script(self,printargs=False,savelog=True):
        """
        Run the external script with the current parameters
        """
        for ps in self.list_ps:
            print("Computing RDF for {} {} ...".format(self.params_input[0],ps.params[self.params_input[0]]))
            ps.run_script(printargs)
        if savelog:
            self.savelog()
        
    def savelog(self):
        """
        Save a log file with all inputs args.
        """
        try:
            os.makedirs("./outputs/logs")
            print("Directory ./outputs/logs Created ")
        except FileExistsError:
            pass
        dateTimeObj = datetime.datetime.now()
        logfilename = dateTimeObj.strftime("%d-%b-%Y_%H-%M-%S.log")
        with open("./outputs/logs/{}".format(logfilename),"w") as logf:
            logf.write("# A line here refers to a line in the input file.\n")
            for ps in self.list_ps:
                dic=vars(ps.args_script)
                line=''.join(['{0} "{1}" '.format(k, v) for k,v in dic.items()])+'\n'
                logf.write(line)                

class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

class ParamSet(object):
    """
    Set of parameters
    """
    def __init__(self,args_script):
        self.compulsory_params=[{"SEL1","REF","TRAJ"},
                                {"SEL1","SIM","FREQ"}]
        self.params={}
        self.args_script=args_script
       # self.args_script=None

    def add_par(self,dic_ext):
        """
        Add parameters with a external dictionary
        """
        for key,value in dic_ext.items():
            self.params[key]=value

    def check_complete(self):
        """
        Check that the whole set of parameters is larger or equal to the list
        of parameters needed for the subscript.
        """
        # Verbosity
        if self.args_script.verbose:
            print(self.params)
            print("Compulsory arguments : {}".format(self.compulsory_params))
        
        current_set=set(self.params.keys())
        
        test,test_missing=[],[]
        for compulsory_set in self.compulsory_params:
            test.append(compulsory_set.issubset(current_set))
            test_missing.append(compulsory_set-current_set)
        if True in test:
            pass
        else :
            assert False,\
            'Missing parameters {} not defined in the config or input file. '.format(test_missing) 
    
    def run_script(self,printargs=True):
        """
        Run external script with attributes given by the current parameters
        """        
        self.args_script=self.get_args()
        if printargs :
            print("Namespace passed to python script: \n {}".format(self.args_script))
        onerdf.run(self.args_script)

    def get_args(self,printargs=False):
        """
        Get the default parameters from the original script
        and replace the values for new parameters defined in ParamSet.params
        We need to provide required=False option for calling runrdf module to 
        prevent help message to stop the script.
        TO DO : doc for format of trajectory files
        """
        # Reset the arguments to avoid interference between args from both python script
        sys.argv=[sys.argv[0]]
        args=onerdf.parse_args(required=False)
        
        # Required arguments
        try :
            args.ref ="{}/{}".format(self.params['PATH'],self.params['REF'])
            args.traj="{}/{}".format(self.params['PATH'],self.params['TRAJ'])
        except :    
            args.ref ="{}/REF_{}.gro".format(self.params['PATH'],self.params['SIM'])
            args.traj="{}/traj_{}_{}pspf.xtc".format(self.params['PATH'],self.params['SIM'],self.params['FREQ'])
        
        args.sel1 =self.params['SEL1'].replace('"','')

        # integer values
        for el in set(self.params.keys()&{'BEGIN','END','NBINS','MAX','YMAX','SUB','FONTSIZE'}):
            setattr(args, el.lower(), int(self.params[el]))

        # boolean values
        for el in set(self.params.keys()&{'CUMULATIVE'}):
            setattr(args, el.lower(), bool(self.params[el]))

        # tuples
        for el in set(self.params.keys()&{'EXCLUSION_BLOCK'}):
            setattr(args, el.lower(), tuple(map(int, re.findall(r'[0-9]+', self.params[el]))))
        
        # string
        for el in set(self.params.keys()&{'SEL2','OUTDIRNAME','SUFFIXNAME'}):
            setattr(args, el.lower(),self.params[el].replace('"',''))
        #if self.args_script.verbose:
        
        if printargs:
            print(args)
        return args

def _format_key(key):
    '''
    Change the format of the parameters to UPPERCASE
    '''
    return key.upper()

# Config File parsing
def read_config(config_file=None,printargs=False):
    """Reads configuration file and returns arguments for `Submitter`.

    Parameters
    ----------
    config_file : str, optional
        Configuration file to use. Default: $HOME/.runRDFrc .

    Returns
    -------
    dict
        Keyword arguments.
    """

    if config_file is None:
        config_file = os.path.expanduser("~/.runRDFrc")
    else:
        config_file = os.path.expanduser(config_file)
    cp = ConfigParser(inline_comment_prefixes=('#', ';'))
    cp.read_file(open(config_file))

    kwargs = {
        "default_parameters": None,
    }    
    # Default Parameters
    if cp.has_section("DefaultParameters"):
        default_parameters = dict()
        for opt in cp.options("DefaultParameters"):
            key = _format_key(opt)
            default_parameters[key] = cp.get("DefaultParameters", opt)
            kwargs["default_parameters"] = default_parameters
    if printargs:
        print(kwargs)
    return kwargs

def parse_args(required=False):
    parser = argparse.ArgumentParser(description="Compute RDF functions recursively on trajectories ,"
                                                "selections, or varying different parameters in"
                                                " default config file and raw input file.",
                                                formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--input_file",
                        default="parameters.in", 
                        help="Inputs file with parameters to apply at each specific trajectory analysis.")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Increase output verbosity")
    parser.add_argument("-c", "--config",
                        default="~/.runRDFrc",
                        help="Config files with default values, for all trajectories analysis.")
    
    args_script = parser.parse_args()
    return parser.parse_args()

def main():
    # parse script arguments
    args_script = parse_args()
    
    # Read defaults parameters
    kwargs=read_config(args_script.config)
    runrdf = RunRdf(args_script,**kwargs)
    
    # Read input file 
    runrdf.read_input(printout=False)

    # Set all parameters for each final script
    runrdf.build_param(**kwargs)
    
    # Check that all parameters needed are provided
    runrdf.check_complete()

    # Run an external script for each set of parameters
    runrdf.run_script(printargs=False,savelog=True)
        
if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description="Compute RDF functions recursively on trajectories ,"
    #                                            "selections, or varying different parameters in"
    #                                            " default config file and raw input file.",
    #                                            formatter_class=RawTextHelpFormatter)
    #parser.add_argument("-i", "--input_file",
    #                    default="parameters.in", 
    #                    help="Inputs file with parameters to apply at each specific trajectory analysis.")
    #parser.add_argument("-v", "--verbose",
    #                    action="store_true",
    #                    help="Increase output verbosity")
    #parser.add_argument("-c", "--config",
    #                    default="~/.runRDFrc",
    #                    help="Config files with default values, for all trajectories analysis.")
    #
    #args_script = parser.parse_args()
    main()