import os
os.system('sleep 5')
xj=590
#seq = [100, 200, 300, 400, 500, 800, 1000, 2000, 4000, 8000]#, 2000, 3000, 4000, 6000, 8000, 10000, 12000, 16000, 20000]
seq = [1,2,4,8,16,32,64]
#alphas = [0.05, 0.01, 0.001, 0.1, 0.2,0.025]
alphas = [0.05]
for jp in seq:
	for a in alphas:
		for x in range(10):
			xj=xj+3
			print './fim_closed mushout'+str(jp)+'-'+str(xj)+'-'+str(a)+' '+str(jp)+' '+str(a)+' mushroom.labels mushroom.trans '+str(xj)
			os.system('./fim_closed mushout'+str(jp)+'-'+str(xj)+'-'+str(a)+' '+str(jp)+' '+str(a)+' mushroom.labels mushroom.trans '+str(xj))
			os.system('sleep 1')