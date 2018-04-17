import luigi
import sciluigi
import os
import subprocess
import logging
import itertools
import pandas as pd
import numpy as np
import ipdb
import matplotlib.pyplot as plt
import seaborn as sns
from puflibs import variables, processing, scltasks

log = logging.getLogger('sciluigi-interface')

class MyWorkflow(sciluigi.WorkflowTask):
    # only required parameter
    outdir = luigi.Parameter()
    cores =  luigi.IntParameter(default=1)
    # genome data
    regions = luigi.Parameter(default='test.bed')

    def workflow(self):
        makebedinput = self.new_task('makebedinput', scltasks.FilenameToTaskOutput, filename=self.regions)

        # divide beds before annotation
        dividebed = self.new_task('dividebed', scltasks.DivideBed10, outdir=os.path.join(self.outdir, 'beds/split'))
        dividebed.in_bed = makebedinput.out_file
        for i in range(2):
            setattr(dividebed, "out_bed%d"%i, sciluigi.TargetInfo(dividebed, dividebed.get_outval(i)))

        

        #return  combinedata
        return dividebed


    
if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow)