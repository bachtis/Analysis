from defs import *

from defs import *

builder = DataSetBuilder(pmap,w,'../../data/JGEN.root','data',1000000000)

data = builder.tree.reduce('massRaw>3.0&&massRaw<3.2')
print 'Estimated mean=',data.mean(data.get().find('massRaw'))
