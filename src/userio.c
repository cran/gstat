/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2011 (C) Edzer Pebesma

    Edzer Pebesma, edzer.pebesma@uni-muenster.de
	Institute for Geoinformatics (ifgi), University of Münster 
	Weseler Straße 253, 48151 Münster, Germany. Phone: +49 251 
	8333081, Fax: +49 251 8339763  http://ifgi.uni-muenster.de 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version. As a special exception, linking 
    this program with the Qt library is permitted.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    (read also the files COPYING and Copyright)
*/

/*
 * userio.c: i/o routines for error, warning, log and progress messages
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "defs.h"

#include "debug.h"
#include "utils.h"
#include "version.h"
#include "userio.h"

void Rprintf(const char *, ...);
void Rf_error(const char *, ...);
#define is_openf(f) (f != NULL)

static FILE *logfile = NULL;

static struct {
	void (*warning_handler)(const char *mess);
	void (*error_handler)(const char *mess, int level);
	void (*printlog_handler)(const char *mess);
	void (*progress_handler)(unsigned int this, unsigned int total);
} gstat_handler = { NULL, NULL, NULL, NULL };

static void (*old_progress_handler)(unsigned int this, unsigned int total) 
		= NULL;

static STRING_BUFFER 
	*error_prefix = NULL, 
	*error_message = NULL, 
	*warning_message = NULL;

static enum Gstat_errno gstat_errno;

const char *error_messages[MAX_ERRNO+1] = {
/* 0 */		"%s",
/* 1 */		"bug in function `%s'",
/* 2 */		"variable not set: %s",
/* 3 */		"variable outside valid range: %s",
/* 4 */		"value not allowed for: %s",
/* 5 */		"no filename set %s",
/* 6 */		"write failed on file `%s'",
/* 7 */		"read failed on file `%s'",
/* 9 */		"cannot read real value from `%s'",
/* 9 */		"cannot read integer from `%s'",
/* 10 */	"syntax error: %s",
/* 11 */	"illegal option or missing argument on `%s'",
/* 12 */	"domain (math) error on `%s'",
/* 13 */	"out of dynamic memory %s",
/* 14 */	"i/o error: %s",
/* 15 */	"no command file%s",
/* 16 */	"%s user interface not compiled in this version",
/* 17 */	"writing to pipe `%s' failed",
/* 18 */	"reading from pipe `%s' failed",
/* 19 */    "function call prevented by secure mode%s",
/* 20 */	"matrix library error: %s",
/* 21 */	"extdbase error: %s"
};

void init_userio(int use_stdio) {
	if (use_stdio) {
		set_gstat_log_file(NULL);
		set_gstat_warning_handler(default_warning);
		set_gstat_error_handler(default_error);
		set_gstat_log_handler(default_printlog);
		set_gstat_progress_handler(default_progress);
	} else {
		/* ... */
	}
	error_prefix    = resize_strbuf(error_prefix, ERROR_BUFFER_SIZE);
	error_message   = resize_strbuf(error_message, ERROR_BUFFER_SIZE);
	warning_message = resize_strbuf(warning_message, ERROR_BUFFER_SIZE);
	error_prefix->str[0] = error_message->str[0] = 
			warning_message->str[0] = '\0';
}

/*
 * error handling function -- print message and error to string, and
 * call error message handler.
 */
void gstat_error(char *fname, int line,
	enum Gstat_errno err_nr, const char *msg) {
	char s[30], *buf;
	int len;

	assert(err_nr <= MAX_ERRNO);
	gstat_errno = err_nr;

	if (error_prefix->str[0] != '\0')
		save_strcat(error_message, error_prefix->str);

	save_strcat(error_message, "gstat: ");
	len = strlen(error_message->str);
	buf = error_message->str + len;
#ifdef HAVE_SNPRINTF
	snprintf(buf, ERROR_BUFFER_SIZE - len,
		error_messages[err_nr], save_string(msg));
#else
	sprintf(buf, error_messages[err_nr], save_string(msg));
#endif
	if (DEBUG_DUMP || err_nr == ER_NULL) { /* print file&line */
		save_strcat(error_message, " (");
		save_strcat(error_message, fname);
		sprintf(s, ", line %d)", line);
		save_strcat(error_message, s);
	}

	if (err_nr == ER_NULL) {
		save_strcat(error_message, "\nVersion info: ");
		/* save_strcat(error_message, GSTAT_OS); */
		save_strcat(error_message, " ");
		save_strcat(error_message, VERSION);
		save_strcat(error_message,
			"\nThis is a bug. Please send the above information, along with\n");
		save_strcat(error_message,
			"the information necessary to reproduce this bug to ");
		save_strcat(error_message, GSTAT_EMAIL);
	}

	gstat_handler.error_handler(error_message->str, err_nr);
	error_message->str[0] = '\0';
	return;
}

/* wrapper function for ErrClo(optopt), in case of error command line option */
void gstat_clo_error(char *f, int l, enum Gstat_errno err, int a) {
	static char s[2];
	sprintf(s, "%c", a);
	gstat_error(f, l, err, s);
} 

/* message() calls for messages preceding a call to ErrMsg() */
void message(char *fmt, ...) {
	va_list args;
	/* char *buf = NULL; */

	va_start(args, fmt);
#ifdef HAVE_VSNPRINTF
	vsnprintf(error_prefix->str, ERROR_BUFFER_SIZE, fmt, args);
#else
	vsprintf(error_prefix->str, fmt, args);
#endif
	va_end(args);
	/* buf = NULL; */
}

/* print a warning message to string, and call warning message handler */
void pr_warning(char *fmt, ...) {
	va_list args;
	char *buf = NULL;

	if (warning_message->max_length < 11)
		resize_strbuf(warning_message, 11);

	warning_message->str[0] = '\0';
	save_strcat(warning_message, "Warning: ");

	buf = warning_message->str + 9;

	va_start(args, fmt);
#ifdef HAVE_VSNPRINTF
	vsnprintf(buf, ERROR_BUFFER_SIZE - 9, fmt, args);
#else
	vsprintf(buf, fmt, args);
#endif
	va_end(args);

	gstat_handler.warning_handler(warning_message->str);
}

void print_progress(unsigned int current, unsigned int total) {
	gstat_handler.progress_handler(current, total);
}

/* get the value of gstat errno */
enum Gstat_errno get_gstat_errno(void) {
	return gstat_errno;
}

/* set the internal gstat errno to NO_ERROR, and reset error mesages */
void reset_gstat_errno(void) {
	assert(error_prefix);
	assert(error_message);

	gstat_errno = ER_NOERROR;
	error_prefix->str[0] = '\0';
	error_message->str[0] = '\0';
}

void set_gstat_warning_handler(void (*warning_fn)(const char *message)) {
	gstat_handler.warning_handler = warning_fn;
}

void set_gstat_error_handler(void (*error_fn)(const char *message, int level)) {
	gstat_handler.error_handler = error_fn;
}

void set_gstat_log_handler(void (*logprint)(const char *str)) {
	gstat_handler.printlog_handler = logprint;
}

void set_gstat_progress_handler(
		void (*progress)(unsigned int this, unsigned int total)) {
	gstat_handler.progress_handler = progress;
}

void push_gstat_progress_handler(
		void (*progress)(unsigned int this, unsigned int total)) {

	assert(old_progress_handler == NULL);

	old_progress_handler = gstat_handler.progress_handler;
	set_gstat_progress_handler(progress);
}

void pop_gstat_progress_handler(void) {

	assert(old_progress_handler != NULL);

	set_gstat_progress_handler(old_progress_handler);
	old_progress_handler = NULL;
}

const char *get_gstat_error_message(void) {
	return (const char *) error_message->str;
}

void print_to_logfile_if_open(const char *mess) {

	if (is_openf(logfile))
		Rprintf("%s", mess);
}

void default_warning(const char *mess) {

	print_to_logfile_if_open(mess);

	Rprintf("%s\n", mess);
	return;
}

void default_error(const char *mess, int level) {

	print_to_logfile_if_open(mess);

	Rf_error("%s\n", mess);
}

void printlog(const char *fmt, ...) {
	STRING_BUFFER *s;
	va_list args;

	s = resize_strbuf(NULL, ERROR_BUFFER_SIZE);

	va_start(args, fmt);
#ifdef HAVE_VSNPRINTF
	vsnprintf(s->str, ERROR_BUFFER_SIZE, fmt, args);
#else
	vsprintf(s->str, fmt, args);
#endif
	va_end(args);

	gstat_handler.printlog_handler(s->str);
	free_strbuf(s);
}

void default_printlog(const char *mess) {

	if (DEBUG_SILENT)
		return;

	if (is_openf(logfile))
		print_to_logfile_if_open(mess);
	else
		Rprintf("%s", mess);
}

int set_gstat_log_file(FILE *f) {
	int retval;

	if (f == NULL) {
		if (is_openf(logfile)) {
			retval = efclose(logfile);
			logfile = NULL;
			return retval;
		} else {
			logfile = NULL;
			return 1;
		}
	} else
		logfile = f;
	return 0;
}

void default_progress(unsigned int current, unsigned int total) {
	static int perc_last = -1, sec_last = -1;
	int perc, sec;
	static time_t start;

	if (total <= 0 || DEBUG_SILENT)
		return;

	if (sec_last == -1) {
		start = time(NULL);
		sec_last = 0;
	}
	perc = floor(100.0 * current / total);
	if (perc != perc_last) { /* another percentage -> calculate time: */
		if (current == total) { /* 100% done, reset: */
			Rprintf("\r%3d%% done\n", 100);
			perc_last = sec_last = -1;
		} else {
			sec = difftime(time(NULL), start);
			if (sec != sec_last) { /* another second -- don't print too often */
				Rprintf("\r%3d%% done", perc);
				perc_last = perc;
				sec_last = sec;
			}
		}
	}
}

void no_progress(unsigned int current, unsigned int total) {
	return;
}
