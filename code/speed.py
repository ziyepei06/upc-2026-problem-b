import numpy as np, time, bkick as B
r0=B.KICK_START; v0=np.array([30.4,9.4,6.9]); om=np.array([0.,0.,-60.])
for dt in (0.004,0.002,0.001,0.0005):
    t=time.time(); fl=B.Flight(cd=0.21)
    tc,y,z,_=fl.goal_plane(r0,v0,om,dt=dt)
    print(dt, "t=%.6f y=%.6f z=%.6f"%(tc,y,z), "dt-wall=%.3fs"%(time.time()-t))
