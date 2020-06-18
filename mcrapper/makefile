all: lcm.c
	gcc -O3 -mavx -o fim_closed -DLCMFREQ_COUNT_OUTPUT lcm.c -lm -w
profile: lcm.c
	gcc -O3 -mavx -o fim_closed -DLCMFREQ_COUNT_OUTPUT -DPROFILE_MINING lcm.c -lm -w
