# SIR
It's a SIR model for novel coronavirus with Eigen which is a C++ library for matrix operations.

You can simulate the situation of the spread of novel coronavirus bys using method createData()in main.cpp. Or you can use your own data to caltulate R0 by using method getRes() which reads data in data document.

I.txt is the daily amount of infected people.
R.txt is the daily amount of recovered people.
S.txt is the amount of susceptible people.

The variable I in method createDate() is the initial number of infected people.
Variable S is the initial number of people who are susceptible.
Variable N is the amount of the people.
Variable R is the initial number of people who are cured.
Variable Pcon is the probability to contact others.
Variable Pspr is the probability to spread the novel coronavirus when people are contated with each other.
Variable Prec is the probability to recover.
Variable iter is number of iterations which means the time of the model.

