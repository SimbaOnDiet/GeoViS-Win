//////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2018, USTC, Liang Yien
//
//  Einstein's Camera - A Relativity Theory Based Physical Renderer
//
//  A Program for the (Undergraduate Course) Optics Presentation
//
//  Modified from GeoVis  
//  And also contained Motion4D, TinyScheme
//  Uses GNU Scientific Library
//  Liang Yien, 
//  Student,
//  University of Science and Technology of China
// 
//  I was confused about the license of Motion4D and GeoVis, in their
//  code they claim they obey the GPL, while the author claim the 
//  software uses CPC License. Anyway, I would be safe to modified
//  the code for presentation in the class. But for the convience
//  of my classmate I might offer them my code , which is done rather
//  privately, like send them code in E-mail.I would immediately stop
//  if I am notified that it violates the author's willing.
//  
//  This code will be kept personal and not distribute publicly given
//  the fact that Motion4D and GeoVis were claimed to applied CPC
//  License. However, I would be gald to have my code in publid freely
//  using GPL if I am convicted that Motion4D and GeoVis is under GPL. 
//
//  Einstein's Camera was  written for my personally purpose only at 
//  the beginning and still is now WITHOUT ANY WARRANTY.No MERCHANTABILITY 
//  or FITNESS FOR A PARTICULAR PURPOSE is guaranteed. 
///////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2013, Universitaet Stuttgart, VISUS, Thomas Mueller
//
//  GeoViS is free software: you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  GeoViS is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with GeoViS.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave
//
//  Motion4D- A library for numerical relativity calulation
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or (at
//  your option) any later version.
//  
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
////////////////////////////////////////////////////////////////////////////////////
// GSL - GNU SCIENTIFIC LIBRARY
// 
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Mark Galassi, James Theiler,
// Brian Gough, Gerard Jungman and many others
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
/////////////////////////////////////////////////////////////////////////////////////
// T I N Y S C H E M E    1 . 3 9
//   Dimitrios Souflis (dsouflis@acm.org)
//   Based on MiniScheme (original credits follow)
// (MINISCM)               coded by Atsushi Moriwaki (11/5/1989)
// (MINISCM)           E-MAIL :  moriwaki@kurims.kurims.kyoto-u.ac.jp
// (MINISCM) This version has been modified by R.C. Secrist.
// (MINISCM)
// (MINISCM) Mini-Scheme is now maintained by Akira KIDA.
// (MINISCM)
// (MINISCM) This is a revised and modified version by Akira KIDA.
// (MINISCM)	current version is 0.85k4 (15 May 1994)
/////////////////////////////////////////////////////////////////////////////////////
// --------------------------------------------------------------------------------
/*
    Mat.h

  Copyright (c) 2009-2014  Thomas Mueller, Frank Grave


   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.

*/

/*!  \class  m4d::Matrix
     \brief  Template class for n x m  matrices.

             n: number of rows, m: number of columns.
*/
// -------------------------------------------------------------------------------- 

#ifndef M4D_MAT_H
#define M4D_MAT_H

#include <cstdio>
#include <iostream>
#include <cstring>
#include <typeinfo>
#include <cassert>

#include "VnD.h"

#ifdef _WIN32
#ifndef __GNUC__
#pragma warning (disable: 4244 )
#endif
#endif

#ifdef _WIN32
#ifdef MATH_EXPORTS
#define MATH_API __declspec(dllexport)
#else /* METRIC_EXPORTS */
#define MATH_API  
#endif /* METRIC_EXPORTS */
#else /* _WIN32 */
#define MATH_API
#endif /* _WIN32 */

namespace m4d
{

//---------------------------------------------------------------------------
//    class-template  mType
//---------------------------------------------------------------------------
template <class mType, int n, int m> class  Matrix
{
protected:
    mType**      mat;
    std::string  classType;
    int          nr,nc;  // number of rows, number of columns
    bool         matIsUnit;

public:
    Matrix ( );
    Matrix ( double val );
    Matrix ( const Matrix &M );
    Matrix ( mType** field );
    ~Matrix();

    void         setAll   ( mType val );
    void         setCoeff ( int row, int col, mType val );
    mType        getCoeff ( int row, int col ) const;

    void         setRow ( int row, const VnD<mType,m> &vec );     // set row-vector
    VnD<mType,m> getRow ( int row ) const;

    void         setCol ( int col, const VnD<mType,n> &vec );     // set col-vector
    VnD<mType,n> getCol ( int col ) const;

    int          getNrRows() const { return nr; }
    int          getNrCols() const { return nc; }

    void         setNull();  // make 0-matrix;
    void         setIdent();  // make identity-matrix

    bool         isIdentMat() const;

    void         transpose();
    void         invert();

    void         getDoubleArray ( double val[] );
    void         getFloatArray  ( float val[] );
    //  mType  trace();

    VnD<mType,m>       operator[]( int z );
    VnD<mType,n>       operator*( const VnD<mType,m> &v );

    void               operator=( const Matrix<mType,n,m> &M );
    Matrix<mType,n,m>  operator+( const Matrix<mType,n,m> &M ) const;
    Matrix<mType,n,m>  operator-( const Matrix<mType,n,m> &M ) const;
    Matrix<mType,n,m>  operator|( const Matrix<mType,n,m> &M ) const;
    Matrix<mType,n,m>  operator*( const mType a ) const;

    Matrix<mType,n,m>  operator^( const int l );

    void               operator*=( const Matrix<mType,n,m> &B );

    void   print  ( std::ostream& ostr = std::cerr ) const;
    void   printS ( FILE* fptr = stderr, const std::string format = "%12.8f " ) const;

    // --------------------- friend functions ---------------------

    // double * Matrix
    friend Matrix<mType,n,m> operator* ( const double a, Matrix<mType,n,m> &M) {
        Matrix<mType,n,m> prod;
        for(int i=0; i<n; i++) {
            for(int j=0; j<m; j++) {
                prod.setCoeff(i,j,a*M.getCoeff(i,j));
                prod.matIsUnit = prod.matIsUnit && (prod.getCoeff(i,j) == double(i==j));
            }
        }
        return prod;
    }

    //  Matrix * VnD (specific matrix vektor multiplication)
    /*
  friend VnD<mType,n> operator* (const Matrix<mType,n,m> &M, const VnD<mType,n> &vec)  {
      VnD<mType,n> q;
      for(int i=0; i<n; i++) {
      q[i] = (mType)0;
      for(int j=0; j<n; j++) {
          q[i] += M.getCoeff(i,j)*vec.x(j);
      }
      }
      return q;
  }
  */
    friend VnD<mType,n> operator* (const Matrix<mType,n,m> &M, const VnD<mType,n> &vec)  {
        VnD<mType,n> q;
        for(int i=0; i<n; i++) {
            q[i] = (mType)0;
            for(int j=0; j<n; j++) {
                q[i] += M.getCoeff(i,j)*vec.x(j);
            }
            q[i]+=M.getCoeff(i,m-1);
        }
        return q;
    }

    // transposeMult  (specific multiplication)
    friend VnD<mType,n> transposeMult ( const Matrix<mType,n,m> &M, const VnD<mType,n> &vec) {
        VnD<mType,n> q;
        for(int i=0; i<n; i++) {
            q[i] = (mType)0;
            for(int j=0; j<n; j++) {
                q[i] += M.getCoeff(j,i)*vec.x(j);
            }
        }
        return q;
    }

    // Matrix * Matrix  ((specific) matrix matrix multiplication)
    friend Matrix<mType,n,m>  operator*(const Matrix<mType,n,m> &m1, const Matrix<mType,n,m> &m2)
    {
        assert(n+1==m || n==m);
        Matrix<mType,n,m> prod;
        mType sum;
        for(int i=0;i<n;i++) {
            for(int j=0;j<n;j++) {
                sum = 0;
                for(int k=0;k<n;k++) {
                    sum += m1.getCoeff(i,k)*m2.getCoeff(k,j);
                }
                prod.setCoeff(i,j,sum);
            }
            // hier geht die spezielle M-M-Multi weiter
            if ((n+1)==m) {
                sum = 0;
                for(int k=0;k<n;k++) {
                    sum += m1.getCoeff(i,k)*m2.getCoeff(k,m-1);
                }
                prod.setCoeff(i,m-1, sum + m1.getCoeff(i,m-1));
            }
        }
        return prod;
    }
};



//---------------------------------------------------------------------------
//      constructor/destructor
//---------------------------------------------------------------------------
template <class mType, int n, int m> Matrix<mType,n,m>::Matrix()
{
    nr=n; nc=m;
    // initialize matrix
    mat = new mType*[n];
    for(int i=0;i<n;i++) {
        mat[i] = new mType[m];
    }
    setNull();  // write null-matrix
}

template <class mType, int n, int m> Matrix<mType,n,m>::Matrix(double val)
{
    nr=n; nc=m;
    // initialize matrix
    mat = new mType*[n];
    for(int i=0;i<n;i++) {
        mat[i] = new mType[m];
    }
    setAll(val);  // write null-matrix
}

template <class mType, int n, int m> Matrix<mType,n,m>::Matrix(mType** field)
{
    nr=n; nc=m;
    // initialize matrix
    mat = new mType*[n];
    for(int i=0;i<n;i++) {
        mat[i] = new mType[m];
    }
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            mat[i][j] = field[i][j];
        }
    }
}

// Copy-constructor
template <class mType, int n, int m> Matrix<mType,n,m>::Matrix(const Matrix<mType,n,m> &M)
{
    mat=new mType*[n];
    for(int i=0;i<n;i++)
        mat[i]=new mType[m];

    int size=sizeof(mType)*m;
    for(int i=0;i<n;i++)
        memcpy(mat[i],M.mat[i],size);
}

template <class mType, int n, int m> Matrix<mType,n,m>::~Matrix()
{
    for(int i=0;i<n;i++) {
        delete [] mat[i];
    }
    delete [] mat;

    mat = NULL;
}


//---------------------------------------------------------------------------
//      set/get coeff()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::setAll(mType val)
{
    for(int i=0;i<n;i++)
        for(int j=0;j<m;j++)
            mat[i][j] = val;
}

template <class mType, int n, int m> void Matrix<mType,n,m>::setCoeff ( int row, int col, mType val )
{
    mat[row][col] = val;
}

template <class mType, int n, int m> mType Matrix<mType,n,m>::getCoeff ( int row, int col ) const
{
    return mat[row][col];
}

//---------------------------------------------------------------------------
//      set/get Row
//---------------------------------------------------------------------------
template <class mType, int n, int m> void  Matrix<mType,n,m>::setRow ( int row, const VnD<mType,m> &vec )
{
    for(int j=0;j<m;j++) {
        mat[row][j] = vec.x(j);
    }
}

template <class mType, int n, int m> VnD<mType,m> Matrix<mType,n,m>::getRow (int row ) const
{
    VnD<mType,m> vec;
    for(int j=0;j<m;j++) {
        vec[j] = mat[row][j];
    }
    return vec;
}

//---------------------------------------------------------------------------
//      set/get Col
//---------------------------------------------------------------------------
template <class mType, int n, int m> void  Matrix<mType,n,m>::setCol ( int col, const VnD<mType,n> &vec )
{
    for(int j=0;j<n;j++) {
        mat[j][col] = vec.x(j);
    }
}

template <class mType, int n, int m> VnD<mType,n> Matrix<mType,n,m>::getCol ( int col ) const
{
    VnD<mType,n> vec;
    for(int j=0;j<n;j++) {
        vec[j] = mat[j][col];
    }
    return vec;
}


//---------------------------------------------------------------------------
//      clear()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::setNull()
{
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            mat[i][j] = (mType)0;
        }
    }
}

//---------------------------------------------------------------------------
//      setIdent()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::setIdent()
{
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            if (i==j) { mat[i][j] = (mType)1; }
            else { mat[i][j] = (mType)0; }
        }
    }
}

//---------------------------------------------------------------------------
//      isIdentMat()
//---------------------------------------------------------------------------
template <class mType, int n, int m> bool Matrix<mType,n,m>::isIdentMat() const
{
    for(int i=0;i<n;i++) {
        for(int j=0;j<n;j++) {
            if ( mat[i][j] != (mType)(i==j) ) return false;
        }
    }
    if (n+1==m) {
        for(int j=0;j<n;j++) {
            if (mat[j][m-1] != (mType)(0)) return false;
        }
    }
    return true;
}

//---------------------------------------------------------------------------
//      transpose()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::transpose()
{
    mType tmp;
    for(int i=0;i<n;i++) {
        for(int j=i;j<m;j++) {
            tmp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = tmp;
        }
    }
}

//---------------------------------------------------------------------------
//      invert()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::invert()
{
    if ((n==m) || (n+1==m)) {
        //std::cerr << "invert..." << std::endl;
        /* Invertierung funktioniert fuer bel. nxn-Matrizen */

        int    i, j, k, r, s, pr, ps, rowj, colj;
        int* row = new int[m];
        int* col = new int[m];
        double Pivot ;

        /* --- Initialisierung --- */

        double** A;
        A = new double*[m];
        for(int i=0;i<m;i++) {
            A[i] = new double[m];
        };

        for ( i = 0 ; i < n ; i++ )
            for ( j = 0 ; j < m ; j++ )
                A[i][j] = mat[i][j];

        if (n+1==m) {
            for(j=0;j<m;j++)
                A[m-1][j] = 0.0;
            A[m-1][m-1] = 1.0;
        }

        /*
    printf("Ausgabe:\n");
    for(i=0;i<m;i++) {
      for(j=0;j<m;j++) {
        printf("%10f ",A[i][j]);
      }
      printf("\n");
    }
    */

        for ( i = 0 ; i < m ; i++ )
            row[i] = col[i] = i ;

        /* --- Austauschverfahren --- */

        for ( j = 0 ; j < m ; j++ )
        {
            /* --- suche Pivot A[row[pr]][col[ps]] --- */

            pr = j ;
            ps = j ;

            for ( r = j ; r < m ; r++ )
                for ( s = j ; s < m ; s++ )
                    if ( fabs(A[row[r]][col[s]]) > fabs(A[row[pr]][col[ps]]) )
                    {
                        pr = r ;   ps = s ;
                    }

            /* --- tausche die Zeilen j,r und die Spalten j,s --- */

            if ( pr > j )
            {
                int temp = row[j];   row[j] = row[pr];   row[pr] = temp;
            }

            if ( ps > j )
            {
                int temp = col[j];   col[j] = col[ps];   col[ps] = temp;
            }

            /* --- im Rest der j-Schleife passiert der eigentliche Austausch
     zwischen der j-ten Zeile und der j-ten Spalte             --- */

            rowj  = row[j] ;
            colj  = col[j] ;
            Pivot = A[rowj][colj] ;

            /* --- 1. Schritt: Rechteckregel --- */

            for ( i = 0 ; i < m ; i++ )
                if ( i != rowj )
                    for ( k = 0 ; k < m ; k++ )
                        if ( k != colj )
                            A[i][k] -= ( A[i][colj] * A[rowj][k] / Pivot ) ;

            /* --- 2. Schritt: Pivotzeile durch negativen Pivot teilen --- */

            for ( k = 0 ; k < m ; k++ )
                if ( k != colj )
                    A[rowj][k] /= (-Pivot) ;

            /* --- 3. Schritt: Pivotspalte durch positiven Pivot teilen --- */

            for ( i = 0 ; i < m ; i++ )
                if ( i != rowj )
                    A[i][colj] /= Pivot ;

            /* --- 4. Schritt: Pivot durch sein Inverses ersetzen --- */

            A[rowj][colj] = 1.0 / Pivot ;

        } /*for j */

        /* inverse Matrix wieder nach mat schreiben */

        for ( i = 0 ; i < m ; i++ )
            for ( k = 0 ; k < m ; k++ )
                if (col[i]<n)
                    mat[col[i]][row[k]] = A[row[i]][col[k]] ;

        for(i=0;i<m;i++) { delete [] A[i];}
        delete [] A;
        delete [] row;
        delete [] col;
    }
    
    // Inverse bei Scale,Translat,Rot
    /*
    if (n+1==m)
      {
    for (int i = 0; i < n; i++ )
      mat[i][3]*=-1.0;
      }
    */
}


//---------------------------------------------------------------------------
//      getDoubleArray / getFloatArray
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::getDoubleArray ( double val[] )
{
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            val[i*m+j] = mat[i][j];
}

template <class mType, int n, int m> void Matrix<mType,n,m>::getFloatArray  ( float val[] )
{
    for(int i=0; i<n; i++)
        for(int j=0; j<m; j++)
            val[i*m+j] = float(mat[i][j]);
}


//---------------------------------------------------------------------------
//      operator[](int i)  //access row vector
//---------------------------------------------------------------------------
template <class mType, int n, int m> VnD<mType,m> Matrix<mType,n,m>::operator[](int z) 
{
    assert(z<nr);
    VnD<mType,m> vec;
    for (int i=0;i<nc;i++) {
        vec.setX(i,mat[z][i]);
    }
    return vec;
}

//---------------------------------------------------------------------------
//      operator=
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::operator=(const Matrix<mType,n,m> &M)
{
    if ( this != &M) {
        for(int i=0;i<nr;i++) {
            delete [] mat[i];
        }
        delete [] mat;
    }

    mat = new mType*[nr];
    for(int i=0;i<nr;i++)
        mat[i] = new mType[nc];

    int size=sizeof(mType)*nc;
    for(int i=0;i<nr;i++) {
        memcpy(mat[i],M.mat[i],size);
    }
}


//---------------------------------------------------------------------------
//      operator+
//---------------------------------------------------------------------------
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator+(const Matrix<mType,n,m> &M) const
{
    Matrix<mType,n,m> Q;
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            Q.mat[i][j] = mat[i][j] + M.mat[i][j];
        }
    }
    return Q;
}

//---------------------------------------------------------------------------
//      operator-
//---------------------------------------------------------------------------
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator-(const Matrix<mType,n,m> &M) const
{
    Matrix<mType,n,m> Q;
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            Q.mat[i][j] = mat[i][j] - M.mat[i][j];
        }
    }
    return Q;
}

//---------------------------------------------------------------------------
//      operator*  
//---------------------------------------------------------------------------
// multiply only matrices of the same type
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator|(const Matrix<mType,n,m> &M) const
{
    assert(n==m);
    Matrix<mType,n,m> Q;

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            Q.mat[i][j] = (mType)0;
            for(int k=0;k<m;k++) {
                Q.mat[i][j] += mat[i][k]*M.mat[k][j];
            }
        }
    }

    return Q;
}


template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator*(const mType a) const
{
    Matrix<mType,n,m> Q;

    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) {
            for(int k=0;k<m;k++) {
                Q.mat[i][j] = a*mat[i][k];
            }
        }
    }
    return Q;
}


//---------------------------------------------------------------------------
//      operator* vector
//---------------------------------------------------------------------------
template <class mType, int n, int m> VnD<mType,n> Matrix<mType,n,m>::operator*(const VnD<mType,m> &v)
{
    VnD<mType,n> q;
    for(int i=0;i<n;i++) {
        q[i] = (mType)0;
        for(int j=0;j<m;j++) {
            q[i] += mat[i][j]*v.x(j);
        }
    }
    return q;
}

//---------------------------------------------------------------------------
//      operator = matrix
//---------------------------------------------------------------------------
/*
template <class mType, int n, int m> const Matrix<mType,n,m>& Matrix<mType,n,m>::operator=( const Matrix<mType,n,m> &B )
{
   if(this != &B)
   {
      int  i ;
      int  j ;
      for(  i = 0; i < 4; i++ )
         for(  j = 0; j < 4; j++ )
            mat[i][j] = B.mat[i][j];

      matIsUnit = B.matIsUnit;
   }

   return *this;
}
*/

//---------------------------------------------------------------------------
//      operator*= matrix
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::operator*=( const Matrix<mType,n,m> &B )
{
    Matrix<mType,n,m> prod( 0.0 );
    prod.matIsUnit = true;

    register int  i ;
    register int  j ;
    register int  k ;
    //int limit;
    
    /*
    // Spezialfall fuer Transformationsmatrizen: Es wird nur die 3x3-Matrix
    // multipliziert.
    if( ( n==3 ) && ( m==4 ))
        limit = 3;
    else
        limit = 4;
        
    for(  i = 0; i < limit; i++ )
        for(  j = 0; j < limit; j++ )
        {
          for(  k = 0; k < limit; k++ )
                prod.mat[i][j] += mat[i][k] * B.mat[k][j];

          prod.matIsUnit = prod.matIsUnit && (prod.mat[i][j] == double(i==j));
        }
    */

    assert( (n==m) || (n+1 == m) );

    if (n==m)
    {
        for( i=0; i<n; i++ )
            for( j=0; j<n; j++ )
            {
                for( k=0; k<n; k++ )
                    prod.mat[i][j] += mat[i][k] * B.mat[k][j];

                prod.matIsUnit = prod.matIsUnit && (prod.mat[i][j] == mType(i==j));
            }
    }
    else if (n+1 == m)
    {
        for( i=0; i<n; i++ )
        {
            for( j=0; j<m; j++ )
            {
                for( k=0; k<n; k++ )
                    prod.mat[i][j] += mat[i][k] * B.mat[k][j];

                prod.matIsUnit = prod.matIsUnit && (prod.mat[i][j] == mType(i==j));
            }
            prod.mat[i][m-1] += mat[i][m-1];
        }
    }

    *this = prod;
}


//---------------------------------------------------------------------------
//      operator^ =matrix
//---------------------------------------------------------------------------
template <class mType, int n, int m> Matrix<mType,n,m> Matrix<mType,n,m>::operator^( const int l )
{
    assert(n==m);
    Matrix<mType,n,m> pow;
    pow.setIdent();

    if (l==0) return pow;

    if (l>0) {
        for(int i=0; i<l; i++) {
            pow = pow*mat;
        }
        pow.matIsUnit = isIdentMat();
        return pow;
    }

    return mat;
}




//---------------------------------------------------------------------------
//      print()
//---------------------------------------------------------------------------
template <class mType, int n, int m> void Matrix<mType,n,m>::print(std::ostream& ostr) const
{
    for(int i=0;i<n;i++) {
        ostr << "( ";
        for(int j=0;j<m; j++) {
            ostr << mat[i][j] << "\t";
        }
        ostr << ")" << std::endl;
    }
}

template <class mType, int n, int m> void  Matrix<mType,n,m>::printS ( FILE* fptr, const std::string format ) const
{  
    for(int i=0;i<n;i++) {
        fprintf(fptr,"(");
        for(int j=0; j<m; j++) {
            fprintf(fptr,format.c_str(),mat[i][j]);
        }
        fprintf(fptr,")\n");
    }
}

} // end namespace m4d

#endif

