import numpy as np, bkick as B
fl=B.FlightBatch(cd=0.21, cl=1.0)
v0=np.tile([30.4,9.4,6.9],(5,1)); om=np.tile([0.,0.,-60.],(5,1))
print(fl.positions_at(v0,om,[0.05,0.10]))
fl2=B.Flight(cd=0.21,cl=1.0)
print(fl2.goal_plane(B.KICK_START,v0[0],om[0]))
print(fl.goal_plane_batch(v0,om))
