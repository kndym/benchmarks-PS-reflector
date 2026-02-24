.PHONY: OBJ

OBJ :
	g++ -O2 \
	  main.cpp \
	  Benchmarks/test_3D_SquareToCircle_logcost_MonteCarlo.h \
	  Benchmarks/test_3D_SquareToTwoGaussSide_logcost_MonteCarlo.h \
	  PushForward/PushForward_Cloud_128.h \
	  PushForward/Pushforward_of_RefRegular.h \
	  QuasiMonteCarlo/Generic_3D_logcost_MonteCarlo.h \
	  QuasiMonteCarlo/MonteCarlo_Pointcloud_3D_128.h \
	  SmallGrid/3D_MonteCarlo_Pointcloud_small.h \
	  SmallGrid/SmallSinkhorn_3D_MC.h \
	  -o main -lm -lpthread

clean :
	$(RM) main


