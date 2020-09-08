Usage:

In python script directory:

`$ python cutout_reproject.py path/masterlist.txt (4400)`

In stack files directory:

`$ python path_to_script/cutout_reproject.py masterlist.txt (4400)`

This script takes masterlist.txt as inlut (which is a list of all the individual stacks from 1ne to 4sw).
It will find the `master_stack_cluster.fits` and make a cutout of it, then reproject all individual stacks to the master cutout.

The second argument 4400 is optional and it is the dimension of master cutout. Default is 4400 pixels.

