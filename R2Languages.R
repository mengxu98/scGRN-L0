

system("/usr/bin/python3.8 TENET ./TENET expression_data.csv 10 trajectory.txt cell_select.txt 1")

system("python TENET ./TENET expression_data.csv 10 trajectory.txt cell_select.txt 1")

system("Rscript SCODE/SCODE.R SCODE/data/exp_train.txt SCODE/data/time_train.txt SCODE/out 100 4 356 100")
