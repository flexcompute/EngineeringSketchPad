/*
 *      CAPS: Computational Aircraft Prototype Syntheses
 *
 *             Phase Utility Application
 *
 *      Copyright 2014-2025, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#define snprintf   _snprintf
#define PATH_MAX   _MAX_PATH
#define SLASH     '\\'
#else
#include <strings.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#define SLASH     '/'
#endif

/* splint fixes */
#ifdef S_SPLINT_S
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
extern ssize_t getline(char ** restrict linep, size_t * restrict linecapp,
                       FILE * restrict stream);
#endif

#include "capsTypes.h"

typedef struct {
  char     *name;               /* the filename */
  CAPSLONG creatim;             /* the file creation time/date stamp */
} dirFiles;


extern /*@null@*/ /*@out@*/ /*@only@*/
       void *EG_alloc(size_t nbytes);
extern /*@null@*/ /*@only@*/
       void *EG_reall(/*@null@*/ /*@only@*/ /*@returned@*/ void *ptr,
                      size_t nbytes);
extern void  EG_free(/*@null@*/ /*@only@*/ void *pointer);
extern /*@null@*/ /*@only@*/
       char *EG_strdup(/*@null@*/ const char *str);

extern int   caps_statFile(const char *path);
extern int   caps_rmFile(const char *path);
extern int   caps_rmDir(const char *path);
extern int   caps_cpDir(const char *src, const char *dst);



#ifdef WIN32
static int getline(char **linep, size_t *linecapp, FILE *stream)
{
  int  i = 0, cnt = 0;

  char *buffer;

  buffer = (char *) EG_reall(*linep, 1026*sizeof(char));
  if (buffer == NULL) return EGADS_MALLOC;

  *linep    = buffer;
  *linecapp = 1024;
  while (cnt < 1024 && !ferror(stream)) {
    i = fread(&buffer[cnt], sizeof(char), 1, stream);
    if (i == 0) break;
    cnt++;
    if ((buffer[cnt-1] == 10) || (buffer[cnt-1] == 13) ||
        (buffer[cnt-1] == 0)) break;
  }

  buffer[cnt] = 0;
  return cnt == 0 ? -1 : cnt;
}
#endif


static void insertLast(int nFile, dirFiles *files)
{
  int      i, j;
  char     *sav;
  CAPSLONG lsav;
  
  if (nFile == 1) return;
  
  lsav = files[nFile-1].creatim;
  sav  = files[nFile-1].name;
  
  for (i = 0; i < nFile-1; i++)
    if (strcmp(files[i].name, sav) > 0) break;
  if (i == nFile-1) return;
  
  for (j = nFile-1; j > i; j--) files[j] = files[j-1];
  files[i].creatim = lsav;
  files[i].name    = sav;
}


static int caps_lsDir(const char *path, int *nFile, dirFiles **files)
{
  int             i;
  short           datetime[6];
  char            full[PATH_MAX];
  CAPSLONG        datim;
  dirFiles        *dirents;
#ifdef WIN32
  WIN32_FIND_DATA ffd;
  HANDLE          hFind;
  SYSTEMTIME      sysTime;
#else
  struct dirent   *de;
  DIR             *dr;
  struct stat     filestat;
  struct tm       tim, *dum;
#endif

  *nFile = 0;
  *files = NULL;
  if (path == NULL) {
    printf(" Information: caps_lsDir called with NULL name!\n");
    return CAPS_NULLNAME;
  }
  
#ifdef WIN32
  snprintf(full, PATH_MAX, "%s\\*", path);
  hFind = FindFirstFile(full, &ffd);
  if (hFind == INVALID_HANDLE_VALUE) {
    printf(" Information: caps_lsDir could not open %s\n", path);
    return CAPS_BADNAME;
  }
  do {
    if (strcmp(ffd.cFileName, ".")  == 0) continue;
    if (strcmp(ffd.cFileName, "..") == 0) continue;
    *nFile += 1;
  } while (FindNextFile(hFind, &ffd) != 0);
  FindClose(hFind);
  if (*nFile == 0) return CAPS_SUCCESS;
  
  *files = dirents = EG_alloc(*nFile*sizeof(dirFiles));
  if (dirents == NULL) return EGADS_MALLOC;
  *nFile = 0;
  hFind  = FindFirstFile(full, &ffd);
  do {
    if (strcmp(ffd.cFileName, ".")  == 0) continue;
    if (strcmp(ffd.cFileName, "..") == 0) continue;
    FileTimeToSystemTime(&ffd.ftCreationTime, &sysTime);
    datetime[0] = sysTime.wYear;
    datetime[1] = sysTime.wMonth;
    datetime[2] = sysTime.wDay;
    datetime[3] = sysTime.wHour;
    datetime[4] = sysTime.wMinute;
    datetime[5] = sysTime.wSecond;
    datim       = datetime[0] - 2000;
    for (i = 1; i < 6; i++) datim = datim*100 + datetime[i];
    dirents[*nFile].name    = EG_strdup(ffd.cFileName);
    dirents[*nFile].creatim = datim;
    *nFile += 1;
    insertLast(*nFile, dirents);
  } while (FindNextFile(hFind, &ffd) != 0);
  FindClose(hFind);
  
#else

  dr = opendir(path);
  if (dr == NULL) {
    printf(" Information: caps_lsDir could not open %s\n", path);
    return CAPS_BADNAME;
  }
  while ((de = readdir(dr)) != NULL) {
    if (strcmp(de->d_name, ".")  == 0) continue;
    if (strcmp(de->d_name, "..") == 0) continue;
    *nFile += 1;
  }
  closedir(dr);
  if (*nFile == 0) return CAPS_SUCCESS;
  
  *files = dirents = EG_alloc(*nFile*sizeof(dirFiles));
  if (dirents == NULL) return EGADS_MALLOC;
  *nFile = 0;
  dr     = opendir(path);
  while ((de = readdir(dr)) != NULL) {
    if (strcmp(de->d_name, ".")  == 0) continue;
    if (strcmp(de->d_name, "..") == 0) continue;
    snprintf(full, PATH_MAX, "%s/%s", path, de->d_name);
    stat(full, &filestat);
/*@-unrecog@*/
    dum = localtime_r(&filestat.st_ctime, &tim); /* this is not the creation */
/*@+unrecog@*/
    if (dum == NULL) {
      datetime[0] = 1900;
      datetime[1] = 0;
      datetime[2] = 0;
      datetime[3] = 0;
      datetime[4] = 0;
      datetime[5] = 0;
    } else {
      datetime[0] = tim.tm_year + 1900;
      datetime[1] = tim.tm_mon  + 1;
      datetime[2] = tim.tm_mday;
      datetime[3] = tim.tm_hour;
      datetime[4] = tim.tm_min;
      datetime[5] = tim.tm_sec;
    }
    datim = datetime[0] - 2000;
    for (i = 1; i < 6; i++) datim = datim*100 + datetime[i];
    dirents[*nFile].name    = EG_strdup(de->d_name);
    dirents[*nFile].creatim = datim;
    *nFile += 1;
    insertLast(*nFile, dirents);
  }
  closedir(dr);
#endif
  
  return CAPS_SUCCESS;
}


static int phaseState(const char *prPath, const char *fileName, int *bFlag)
{
  int  stat;
  char full[PATH_MAX];
  
  *bFlag = 0;
  snprintf(full, PATH_MAX, "%s%c%s", prPath, SLASH, fileName);
  stat = caps_statFile(full);
  if (stat != EGADS_OUTSIDE) return EGADS_NOTFOUND;
  snprintf(full, PATH_MAX, "%s%c%s%ccapsRestart",
           prPath, SLASH, fileName, SLASH);
  stat = caps_statFile(full);
  if (stat != EGADS_OUTSIDE) return EGADS_NOTFOUND;
  
  snprintf(full, PATH_MAX, "%s%c%s%ccapsLock",
           prPath, SLASH, fileName, SLASH);
  stat = caps_statFile(full);
  if (stat == EGADS_SUCCESS) *bFlag += 1;

  snprintf(full, PATH_MAX, "%s%c%s%ccapsClosed",
           prPath, SLASH, fileName, SLASH);
  stat = caps_statFile(full);
  if (stat == EGADS_SUCCESS) *bFlag += 2;
  
  return EGADS_SUCCESS;
}


static void listPhases(const char *prPath)
{
  int      i, j, n, stat, nFile, flag;
  char     full[PATH_MAX];
  dirFiles *files;
  FILE     *fp;
  
  stat = caps_lsDir(prPath, &nFile, &files);
  if (stat != CAPS_SUCCESS) return;
  
  printf("                   Phase Name:                    Parent Name:\n");
  printf("                   -------------------------------------------\n");
  for (i = 0; i < nFile; i++) {
    stat = phaseState(prPath, files[i].name, &flag);
    if (stat != EGADS_SUCCESS) {
      EG_free(files[i].name);
      continue;
    }
    if ((flag&1) == 1) {
      printf(" Locked  ");
    } else {
      printf("         ");
    }
    if ((flag&2) == 2) {
      printf(" Closed  ");
    } else {
      printf("         ");
    }
    n = printf(" %s", files[i].name);
    snprintf(full, PATH_MAX, "%s%c%s%cparent.txt", prPath, SLASH, files[i].name,
             SLASH);
    EG_free(files[i].name);
    fp = fopen(full, "r");
    if (fp == NULL) {
      printf("\n");
    } else {
      fscanf(fp, "%s", full);
      fclose(fp);
      if (full[0] == '0') {
        printf("\n");
      } else {
        for (j = n; j < 32; j++) printf(" ");
        printf("%s\n", full);
      }
    }
  }
  EG_free(files);
  printf("\n");
}


static void showLock(const char *prPath, const char *phName)
{
  int    stat, flag;
  char   full[PATH_MAX];
  size_t linecap = 0;
  char   *line = NULL;
  FILE   *fp;

  stat = phaseState(prPath, phName, &flag);
  if (stat == EGADS_NOTFOUND) {
    printf(" %s%c%s is not a valid Phase!\n", prPath, SLASH, phName);
    return;
  }
  if ((flag&1) != 1) {
    printf(" Phase %s%c%s is not Locked!\n", prPath, SLASH, phName);
    return;
  }
  
  snprintf(full, PATH_MAX, "%s%c%s%ccapsLock", prPath, SLASH, phName, SLASH);
  fp = fopen(full, "r");
  if (fp == NULL) {
    printf(" Cannot open Lock file: %s!\n", full);
    return;
  }
  (void) getline(&line, &linecap, fp);
  fclose(fp);
  printf(" Lock File => %s\n", line);
  EG_free(line);
}


static void removeLock(const char *prPath, const char *phName)
{
  int  stat, flag;
  char full[PATH_MAX];

  stat = phaseState(prPath, phName, &flag);
  if (stat == EGADS_NOTFOUND) {
    printf(" %s%c%s is not a valid Phase!\n", prPath, SLASH, phName);
    return;
  }
  if ((flag&1) != 1) {
    printf(" Phase %s%c%s is not Locked!\n", prPath, SLASH, phName);
    return;
  }
  
  snprintf(full, PATH_MAX, "%s%c%s%ccapsLock", prPath, SLASH, phName, SLASH);
  stat = caps_rmFile(full);
  if (stat != EGADS_SUCCESS) {
    printf(" Cannot remove %s = %d\n", full, stat);
  } else {
    printf(" Lock File removed!\n\n");
  }
}


static void deletePhase(const char *prPath, const char *phName)
{
  int      i, j, k, len, stat, nFile, nPFile, flag, phlen;
  char     full[PATH_MAX], lnk[PATH_MAX], src[PATH_MAX];
  dirFiles *files, *pfiles;
  FILE     *fp;

  stat = phaseState(prPath, phName, &flag);
  if (stat == EGADS_NOTFOUND) {
    printf(" %s%c%s is not a valid Phase!\n", prPath, SLASH, phName);
    return;
  }
  phlen = strlen(phName);

  /* look at all other Phases */
  stat = caps_lsDir(prPath, &nFile, &files);
  if (stat != CAPS_SUCCESS) return;
  
  /* do we have parents that point to us? */
  for (k = i = 0; i < nFile; i++) {
    if (strcmp(files[i].name, phName) == 0) continue;
    
    stat = phaseState(prPath, files[i].name, &flag);
    if (stat != EGADS_SUCCESS) continue;
    
    snprintf(full, PATH_MAX, "%s%c%s%cparent.txt",
             prPath, SLASH, files[i].name, SLASH);
    fp = fopen(full, "r");
    if (fp == NULL) continue;
    fscanf(fp, "%s", lnk);
    fclose(fp);
    if (strcmp(lnk, phName) != 0) continue;
    k++;
    break;
  }
  if (k != 0) {
    printf(" %s NOT deleted -- parent of %s!\n", phName, files[i].name);
    goto cleanup;
  }

  /* find any with links that may point to us */
  for (i = 0; i < nFile; i++) {
    if (strcmp(files[i].name, phName) == 0) continue;

    stat = phaseState(prPath, files[i].name, &flag);
    if (stat != EGADS_SUCCESS) continue;

    snprintf(full, PATH_MAX, "%s%c%s", prPath, SLASH, files[i].name);
    stat = caps_lsDir(full, &nPFile, &pfiles);
    if (stat != CAPS_SUCCESS) continue;
    
    /* look for clnk files*/
    for (j = 0; j < nPFile; j++) {
      if (pfiles[j].name == NULL) continue;
      len = strlen(pfiles[j].name);
      if (len > 5)
        if ((pfiles[j].name[len-5] == '.') && (pfiles[j].name[len-4] == 'c') &&
            (pfiles[j].name[len-3] == 'l') && (pfiles[j].name[len-2] == 'n') &&
            (pfiles[j].name[len-1] == 'k')) {
          snprintf(full, PATH_MAX, "%s%c%s%c%s",
                   prPath, SLASH, files[i].name, SLASH, pfiles[j].name);
          fp = fopen(full, "r");
          if (fp == NULL) {
            printf(" Cannot open %s!\n", full);
          } else {
            fscanf(fp, "%s", lnk);
            fclose(fp);
            if (strlen(lnk) > phlen) {
              for (k = 0; k < phlen; k++)
                if (lnk[k] != phName[k]) break;
              if ((k == phlen) && (lnk[phlen] == SLASH)) {
                /* found a link to us! */
                len = strlen(full);
                full[len-5] = 0;
                snprintf(src, PATH_MAX, "%s%c%s", prPath, SLASH, lnk);
                printf(" %s repopulated from %s\n", full, src);
                stat = caps_cpDir(src, full);
                if (stat != EGADS_SUCCESS)
                  printf(" Cannot copy directory = %d\n", stat);
                full[len-5] = '.';
                stat = caps_rmFile(full);
                if (stat != EGADS_SUCCESS)
                  printf(" Cannot remove %s = %d\n", full, stat);
              }
            }
          }
        }
      EG_free(pfiles[j].name);
    }
    EG_free(pfiles);
  }
  
  snprintf(full, PATH_MAX, "%s%c%s", prPath, SLASH, phName);
  stat = caps_rmDir(full);
  if (stat != EGADS_SUCCESS) {
    printf(" Cannot remove %s = %d\n", full, stat);
  } else {
    printf(" %s deleted!\n", full);
  }

cleanup:
  for (i = 0; i < nFile; i++) EG_free(files[i].name);
  EG_free(files);
  printf("\n");
}


static void fullCopy(const char *prPath, const char *newPath)
{
  int stat;

  stat = caps_statFile(newPath);
  if (stat != EGADS_NOTFOUND) {
    printf(" Destination %s exists!\n", newPath);
    return;
  }
  stat = caps_cpDir(prPath, newPath);
  if (stat != EGADS_SUCCESS) {
    printf(" Cannot do a full copy of %s to %s = %d\n",
           prPath, newPath, stat);
  } else {
    printf(" Full copy of %s to %s complete!\n\n", prPath, newPath);
  }
}


static void copyPhase(const char *prPath, const char *phName,
                      const char *newPath)
{
  int      i, len, stat, flag, nFile;
  char     dst[PATH_MAX], src[PATH_MAX], full[PATH_MAX], lnk[PATH_MAX];
  dirFiles *files;
  FILE     *fp;
  
  stat = phaseState(prPath, phName, &flag);
  if (stat == EGADS_NOTFOUND) {
    printf(" %s%c%s is not a valid Phase!\n", prPath, SLASH, phName);
    return;
  }
  snprintf(src, PATH_MAX, "%s%c%s", prPath,  SLASH, phName);
  snprintf(dst, PATH_MAX, "%s%c%s", newPath, SLASH, phName);
  stat = caps_statFile(dst);
  if (stat != EGADS_NOTFOUND) {
    printf(" Destination %s exists!\n", dst);
    return;
  }
  
  stat = caps_cpDir(src, dst);
  if (stat != EGADS_SUCCESS) {
    printf(" Cannot do a copy of %s to %s = %d\n", src, dst, stat);
    return;
  }
  /* remove parent info if exists */
  snprintf(lnk, PATH_MAX, "%s%cparent.txt", dst, SLASH);
  (void) caps_rmFile(lnk);
  
  /* adjust any links */
  stat = caps_lsDir(dst, &nFile, &files);
  if (stat != CAPS_SUCCESS) return;
  for (i = 0; i < nFile; i++) {
    if (files[i].name == NULL) continue;
    len = strlen(files[i].name);
    if (len > 5)
      if ((files[i].name[len-5] == '.') && (files[i].name[len-4] == 'c') &&
          (files[i].name[len-3] == 'l') && (files[i].name[len-2] == 'n') &&
          (files[i].name[len-1] == 'k')) {
        /* remove the link */
        snprintf(full, PATH_MAX, "%s%c%s", dst, SLASH, files[i].name);
        fp = fopen(full, "r");
        if (fp == NULL) {
          printf(" Cannot open %s!\n", full);
        } else {
          fscanf(fp, "%s", lnk);
          fclose(fp);
          snprintf(src, PATH_MAX, "%s%c%s", prPath, SLASH, lnk);
          len = strlen(full);
          full[len-5] = 0;
          printf(" %s repopulated from %s\n", full, src);
          stat = caps_cpDir(src, full);
          if (stat != EGADS_SUCCESS)
            printf(" Cannot copy directory = %d\n", stat);
          full[len-5] = '.';
          stat = caps_rmFile(full);
          if (stat != EGADS_SUCCESS)
            printf(" Cannot remove %s = %d\n", full, stat);
        }
      }
    EG_free(files[i].name);
  }
  EG_free(files);
  
  printf(" Copy of %s to %s complete!\n\n", src, dst);
}


static void phaseProblem(const char *prPath, const char *phName,
                         const char *newName)
{
  int  stat, flag;
  char dst[PATH_MAX], src[PATH_MAX];
  
  stat = phaseState(prPath, phName, &flag);
  if (stat == EGADS_NOTFOUND) {
    printf(" %s%c%s is not a valid Phase!\n", prPath, SLASH, phName);
    return;
  }
  snprintf(src, PATH_MAX, "%s%c%s", prPath, SLASH, phName);
  snprintf(dst, PATH_MAX, "%s%c%s", prPath, SLASH, newName);
  stat = caps_statFile(dst);
  if (stat != EGADS_NOTFOUND) {
    printf(" Destination %s exists!\n", dst);
    return;
  }
  
  stat = caps_cpDir(src, dst);
  if (stat != EGADS_SUCCESS) {
    printf(" Cannot do a copy of %s to %s = %d\n", src, dst, stat);
  } else {
    printf(" Copy of %s to %s complete!\n\n", src, dst);
  }
}
  

int main(int argc, char *argv[])
{
  int hit = 0;
  
  printf(" phaseUtil for ESP Rev %d.%02d\n\n", CAPSMAJOR, CAPSMINOR);
  
  if (argc == 2) {
    if (argv[1][0] != '-') {
      listPhases(argv[1]);
      hit++;
    }
  } else if (argc == 4) {
    if (argv[1][0] != '-') {
      if (strcmp(argv[2], "-l") == 0) {
        showLock(argv[1], argv[3]);
        hit++;
      } else if (strcmp(argv[2], "-r") == 0) {
        removeLock(argv[1], argv[3]);
        hit++;
      } else if (strcmp(argv[2], "-d") == 0) {
        deletePhase(argv[1], argv[3]);
        hit++;
      } else if (strcmp(argv[2], "-f") == 0) {
        fullCopy(argv[1], argv[3]);
        hit++;
      }
    }
  } else if (argc == 5) {
    if (argv[1][0] != '-') {
      if (strcmp(argv[2], "-c") == 0) {
        copyPhase(argv[1], argv[3], argv[4]);
        hit++;
      } else if (strcmp(argv[2], "-p") == 0) {
        phaseProblem(argv[1], argv[3], argv[4]);
        hit++;
      }
    }
  }
  
  if (hit == 0) {
    printf(" Usage: phaseUtil prPath                   -- list Phases\n");
    printf("        phaseUtil prPath -l phName         -- show lock owner\n");
    printf("        phaseUtil prPath -r phName         -- remove lock\n");
    printf("        phaseUtil prPath -c phName newPath -- copy Phase\n");
    printf("        phaseUtil prPath -p phName newName -- copy Phase in Problem\n");
    printf("        phaseUtil prPath -d phName         -- delete Phase\n");
    printf("        phaseUtil prPath -f newPath        -- full copy of Problem");
    printf("\n\n");
    return 1;
  }
  
  return 0;
}
