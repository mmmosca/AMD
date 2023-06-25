/*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*/
#include "cmd.h"

int isCharInString(char c, char* str) {
		int i = 0, len=strlen(str);
		while (i < len) {
			if (str[i] == c) return i;/* return the index */
			i++;
		}
		return -1; /*not found*/
	}

int AreStringsEqualFrom(const char* s1, const char* s2,int from) {
		int len1, len2, i=0, lensub, lenbig;
		const char* bigger, *smaller;

		if ((s1 == NULL) || (s2 == NULL)) return 0;

		len1 = strlen(s1);
		len2 = strlen(s2);

		if (len1 > len2) {
			bigger = s1;
			smaller = s2;
			lensub = len2;
			lenbig = len1;
		}
		else
		{
			bigger = s2;
			smaller = s1;
			lensub = len1;
			lenbig = len2;
		}

		if ((lenbig - from) < lensub) return 0;

		while ((bigger[i+from] == smaller[i]) && (i < lensub)) {
			i++;
		}

		if (i == lensub) return 1;
		else return 0;
	}

int isSubstring(char* sub, char* str) {
		int sublen= strlen(sub), len = strlen(str), i=0;

		while (i <= len-sublen) {
			if (AreStringsEqualFrom(sub, str,0)) {
				return i;
			}
			i++;
			str++;
		}

		return -1;
	}

void CommandLine::reset_values() {
		optind = 1;
		argvind = 2;
		myoptarg = NULL;
		formatind = -1;
		argformatind = 0;
	}

char CommandLine::mygetopt(int argc, char** argv, char* format) {
	
	/*if optind is one means that there are no more options*/
	if (optind == -1) {
		start = 1;		
		return -1;
	}
	/*start from 1 and skip the name of the executable in the command line*/
	if (start) {
		reset_values();
		formatind = -1;
		argformatind = 0;
		start = 0;
	}

	/*if optind goes over the number of argument commandline*/
	if (optind > argc-1) {
		PRINT_DATA();
		puts("Too options!");
		exit(EXIT_FAILURE);
	}
	/*option element must have 2 characters (i.e. -g)*/
	if (strlen(argv[optind]) != 2) {
		PRINT_DATA();
		puts("Options syntax invalid!");
		exit(EXIT_FAILURE);
	}
	/*the first element of an option must be '-'*/
	if (argv[optind][0] != '-') {
		PRINT_DATA();
		puts("It is not an option!");
		exit(EXIT_FAILURE);
	}

	if ((formatind = isCharInString(argv[optind][1], format)) >= 0) {
		argformatind = formatind + 1;
		/*Format char is at the end of the format line*/
		if (argformatind < (int) strlen(format)) {
			argvind = optind + 1;
			/*if the argument index (argvind) in the array is not out of bounds*/
			if (argvind <= argc-1) {
				/*if there is a ':' in format string, argument is expected after*/
				if (format[argformatind] == ':') {
					/*argument expected, update the myoptarg*/
					if (argv[argvind][0] != '-') {
						myoptarg = argv[argvind];
#ifdef DEBUG
						PRINT_DATA();
#endif // DEBUG
						optind += 2;
					}
					/*there is '-' in the argument*/
					else
					{
						PRINT_DATA();
						puts("This option must have an argument!");
						exit(EXIT_FAILURE);
					}
				}
				/*No argument expected for optind option*/
				else {
					/*check forward if there is no option*/
					if (argv[argvind][0] != '-') {
						PRINT_DATA();
						puts("This option must not have an argument!");
						exit(EXIT_FAILURE);
					}
					myoptarg = NULL;
#ifdef DEBUG
					PRINT_DATA();
#endif // DEBUG
					optind++;
				}
			}
			/*if the argument index (argvind) in the array is out of bounds*/
			else {
				myoptarg = NULL;
#ifdef DEBUG
				PRINT_DATA();
#endif // DEBUG
				optind++;
			}
			
		}
		else {
			myoptarg = NULL;
#ifdef DEBUG
			PRINT_DATA();
#endif // DEBUG
			optind++;
		}


		if (optind > argc-1) optind = -1;
		return format[formatind];
	}
	else {
		myoptarg = NULL;
		PRINT_DATA();
		puts("Unknown option!");
		exit(EXIT_FAILURE);
	}	
}

char* CommandLine::mygetoptW(int argc, char** argv, char* format) {

		/*if optind is one means that there are no more options*/
		if (optind == -1) {
			start = 1;
			return NULL;
		}
		/*start from 1 and skip the name of the executable in the command line*/
		if (start) {
			reset_values();
			start = 0;
		}

		/*if optind goes over the number of argument commandline*/
		if (optind > argc - 1) {
			PRINT_DATAW();
			puts("Error in parameters number!");
			return NULL;
		}
		/*option element must have more than 2 characters (i.e. -steps)*/
		if (strlen(argv[optind]) <= 2) {
			PRINT_DATAW();
			puts("Options syntax invalid!");
			exit(EXIT_FAILURE);
		}
		/*the first element of an option must be '-'*/
		if (argv[optind][0] != '-') {
			PRINT_DATAW();
			puts("It is not an option!");
			exit(EXIT_FAILURE);
		}

		curr_option = ++argv[optind];
		if ((formatind = isSubstring(argv[optind], format)) >= 0) {
			argformatind = formatind + strlen(argv[optind]);
			/*Format substring is at the end of the format line*/
			if (argformatind < (int)strlen(format)) {
				argvind = optind + 1;
				/*if the argument index (argvind) in the array is not out of bounds*/
				if (argvind <= argc - 1) {
					/*if there is a ':' in format string, argument is expected after*/
					if (format[argformatind] == ':') {
						/*argument expected, update the myoptarg*/
						if (argv[argvind][0] != '-') {
							myoptarg = argv[argvind];
	#ifdef DEBUG
							PRINT_DATAW();
	#endif // DEBUG
							optind += 2;
						}
						/*there is '-' in the argument*/
						else
						{
							myoptarg = argv[argvind];
							PRINT_DATAW();
							puts("This option must have a parameter!");
							exit(EXIT_FAILURE);
						}
					}
					/*No argument expected for optind option*/
					else {

						/*check forward if there is no option*/
						if (argv[argvind][0] != '-') {
							myoptarg = argv[argvind];
							PRINT_DATAW();
							puts("This option must not have a parameter!");
							exit(EXIT_FAILURE);
						}

						if (format[argformatind] != '|') {
							PRINT_DATAW();
							puts("Format and input are not valid!");
							exit(EXIT_FAILURE);
						}
						myoptarg = NULL;
	#ifdef DEBUG
						PRINT_DATAW();
	#endif // DEBUG
						optind++;
					}
				}
				/*if the argument index (argvind) in the array is out of bounds*/
				else {
					/*No arguments expected in format*/
					if (format[argformatind] != '|') {
						PRINT_DATAW();
						puts("Format and input are not valid!");
						exit(EXIT_FAILURE);
					}
					myoptarg = NULL;
	#ifdef DEBUG
					PRINT_DATAW();
	#endif // DEBUG
					optind++;
				}

			}
			/*Format substring is at the end of the format line*/
			else {
				if (format[argformatind] != '|') {
					PRINT_DATAW();
					puts("Format and input are not valid!");
					exit(EXIT_FAILURE);
				}
				myoptarg = NULL;
	#ifdef DEBUG
				PRINT_DATAW();
	#endif // DEBUG
				optind++;
			}

			if (optind > argc - 1) optind = -1;
			return curr_option;
		}
		else {
			myoptarg = NULL;
			PRINT_DATAW();
			puts("Unknown option!");
			exit(EXIT_FAILURE);
		}
	}

char* strsep(char** elem_pointer, char* pattern) {
		char* current_tkn, *p;
		int i = 0, j, len;

		if (*elem_pointer == NULL || pattern == NULL) exit(EXIT_FAILURE);	
		len = strlen(*elem_pointer); /*len of the entire string*/
		p = *elem_pointer;
		while ((isCharInString(p[i], pattern) == -1) && (i < len)) {
			i++;
		}
		len = i; /*len of the first token*/
		current_tkn = (char*) malloc(len * sizeof(char)+1);
	
		j = 0;
		for (i = 0; i < len; i++) {
			if (p[i] != ' ') {
				current_tkn[j] = p[i];
				j++;
			}	
		}
		current_tkn[j] = '\0';

		*elem_pointer += (len+1);
		return current_tkn;
	}
