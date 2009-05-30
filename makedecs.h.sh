#!/bin/sh
sed 's/SFTYPE /extern SFTYPE /' defs.h > decstemp1.h
sed 's/FTYPE /extern FTYPE /' decstemp1.h > decstemp2.h
sed 's/extern Sextern FTYPE /extern SFTYPE /' decstemp2.h > decstemp3.h
sed 's/int /extern int /' decstemp3.h > decstemp4.h
sed 's/short /extern short /' decstemp4.h > decstemp5.h
sed 's/MPI_Status /extern MPI_Status /' decstemp5.h > decstemp6.h
sed 's/MPI_Group /extern MPI_Group /' decstemp6.h > decstemp7.h
sed 's/MPI_Comm /extern MPI_Comm /' decstemp7.h > decstemp8.h
sed 's/char /extern char /' decstemp8.h > decstemp9.h
sed 's/FILE\*/extern FILE\* /' decstemp9.h > decstemp10.h
sed 's/mpidefs.h/mpidecs.h/' decstemp10.h > decstemp11.h
sed 's/double /extern double /' decstemp11.h > decstemp12.h
sed 's/float /extern float /' decstemp12.h > decstemp13.h
sed 's/long /extern long /' decstemp13.h > decs.h
rm decstemp*.h
