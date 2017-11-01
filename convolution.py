import cmath
import numpy as np
from itertools import islice , chain
from random import uniform

def compMult( v1 , v2 ):
    result = []
    for j in range(len(v1)):
        result += [ (v1[j])*(v2[j]) ]
    return result

def polyEval( A, w ):
    result = 0
    for j in range(len(A)):
        result += (A[j])*(w**j)
    return result

def getOmega( n ):
    return np.exp( (2*np.pi*np.complex(0,1)) / n )

def FFT( A ):
    w = getOmega(len(A))
    return FFTlocal( A , w )

def FFTlocal( A , w ):
    if len(A) == 1:
        return A
    n = len(A)
    Ae = FFTlocal( A[0::2] , w**2 )
    Ao = FFTlocal( A[1::2] , w**2 )
    r = [np.complex(0,0)] * len(A)
    for j in range( int((n/2)) ):
        r[j] = ((Ae[j] + (w**j)*Ao[j]))
        r[j + int(n/2)] = ((Ae[j] - (w**j)*Ao[j]))
    return r

#[round(this,5) for this in r]

def invFFT( x ):
    reals = [np.real(a)*(1/4) for a in x]
    reals1 = reals[::-1]
    realcar = reals1.pop()
    reals1.insert(0,realcar)
    imags = [np.imag(a)*(1/4) for a in x]
    imags1 = imags[::-1]
    imagscar = imags.pop()
    imags1.insert(0,imagscar)
    return FFT([np.complex(r,i) for (r,i) in zip(reals1, imags1)])

def invSFT( x ):
    reals = [np.real(a)*(1/4) for a in x]
    reals1 = reals[::-1]
    realcar = reals1.pop()
    reals1.insert(0,realcar)
    imags = [np.imag(a)*(1/4) for a in x]
    imags1 = imags[::-1]
    imagscar = imags.pop()
    imags1.insert(0,imagscar)
    return SFT([np.complex(r,i) for (r,i) in zip(reals1, imags1)])
    
def SFT( A ):
    n = len(A)
    w = getOmega( n )
    res = []
    for j in range( n ):
        res.append( round( polyEval( A , w**j ) , 5 ) )
    return res

def convolutionFast( V1 , V2 ):
    for j in range( len(V1)):
        V1 += [0]
        V2 += [0]
    FFTv1 = FFT( V1)
    FFTv2 = FFT( V2)
    compMultRes = compMult( FFTv1, FFTv2)
    cmplxRes = invFFT( compMultRes )
    return [ np.real( this) for this in cmplxRes ]

def convolutionSlow( V1 , V2 ):
    for j in range( len(V1)):
        V1 += [0]
        V2 += [0]
    SFTv1 = SFT( V1)
    SFTv2 = SFT( V2)
    compMultRes = compMult( SFTv1, SFTv2)
    cmplxRes = invSFT( compMultRes )
    return [ np.real( this) for this in cmplxRes ]

def genVectors( n ):
    v1 , v2 = [] , []
    for j in range(n):
        v1 += [ round( uniform( 0, 1), 1) ]
        v2 += [ round( uniform( 0, 1), 1) ]
    return ( v1, v2 )








def main():
    ( vec1, vec2) = genVectors(5)
    res1 = convolutionSlow( vec1, vec2)
    print( res1 ) 
    '''
    moreInput = True
    while moreInput:
        print( 'Enter integer n for size of vector:')
        n = input('you have sleceted:')
        ( vec1, vec2) = genVectors( n)
        if n <= 100:
            resFast = convolutionFast( vec1, vec2)
            resSlow = convolutionSlow( vec1, vec2)
    '''

if __name__ == "__main__":

    main()