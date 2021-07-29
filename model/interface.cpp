#include "header.h"


void interface::getRunOption(string addressRUN, string addressOPT)
{
 	int a;
	cout << "for a single run press 0, for optimization press 1:";
	cin >> runOption;
	while(runOption != 0 && runOption != 1)
	{
		cout << endl<< " 0 or 1?";
		cin >> runOption;
	}
	if(runOption == 0)
	{
		cout << endl << "I assume that you have set the parameters (D, K, m, n, and Tc) and the initial file name in RUN.C";
		cout << endl <<endl<< "the current output folder is:" << addressRUN;
	}
	if(runOption == 1)
	{
		cout << endl << " I assume that you have set optimization parametersp (pop and iter) in optimization.C";
		cout << endl<< "the current output folder is:" << addressOPT;
	}
	cout << endl<< "Press any integer number to continue ;)";
	cin >> a;
}
