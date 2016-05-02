import sys
import structopt

input_filename = sys.argv[1]
structopt.setup(input_filename)
optimizer = structopt.Optimizer()
optimizer.run()
