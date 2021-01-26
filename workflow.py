from gwf import Workflow, AnonymousTarget
import os
import glob
import pandas as pd
import genominterv
from groups import Group

gwf = Workflow(defaults={"account": "baboondiversity"})

