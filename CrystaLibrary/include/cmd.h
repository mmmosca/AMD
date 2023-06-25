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
#ifndef CMD_H_
#define CMD_H_

#ifdef DEBUG
#define PRINT_DATA()	printf("CURRENT ---- OPTIND: %d ---- ARGVIND: %d ---- OPTION: %c ---- OPTARG: %s\n", optind, argvind, argv[optind][1], myoptarg);
#define PRINT_DATAW()	printf("CURRENT ---- OPTIND: %d ---- ARGVIND: %d ---- OPTION: %s ---- OPTARG: %s\n", optind, argvind, argv[optind], myoptarg);
#else
#define PRINT_DATA()	printf("OPTION: %c ---- OPTARG: %s\n", argv[optind][1], myoptarg);
#define PRINT_DATAW()	printf("OPTION: %s ---- OPTARG: %s\n", argv[optind], myoptarg);
#endif

#include <stdio.h>
#include <string.h>
#include <Windows.h>

int isCharInString(char c, char* str);

int AreStringsEqualFrom(const char* s1, const char* s2, int from);

int isSubstring(char* sub, char* str);

char* strsep(char** elem_pointer, char* pattern);

struct CommandLine {
private:
	/*option index in the array of line commands*/
	/*argument index in the array of line commands*/
	int optind, argvind, start, formatind, argformatind;
	/*pointer to the argument*/
	char *curr_option;
public:
	char* myoptarg;

	CommandLine() : start{ 1 }, formatind{ -1 }, argformatind{ 0 } {
		reset_values();
	};

	void reset_values();

	char mygetopt(int argc, char** argv, char* format);

	char* mygetoptW(int argc, char** argv, char* format);
};
#endif /* CMD_H_ */