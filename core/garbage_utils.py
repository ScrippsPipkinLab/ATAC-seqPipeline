


#### You do not need this file. This is garbage. 

#### Configure controls and treatment 
import itertools
import numpy 

c = numpy.unique((sampsheet[sampsheet['C/T']=='C']['Status']).to_numpy())
t = numpy.unique((sampsheet[sampsheet['C/T']=='T']['Status']).to_numpy())

combs = list(itertools.product(c, t))

for control, treatment in combs: 
    control_files = (sampsheet[sampsheet['Status']==control]['SampleName']).to_numpy()
    treatment_files = (sampsheet[sampsheet['Status']==treatment]['SampleName']).to_numpy()

    print(' '.join(control_files + '.bam'))
    print(' '.join(treatment_files + '.bam'))