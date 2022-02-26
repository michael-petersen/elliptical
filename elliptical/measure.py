# different measures for ellipses


class measureEllipse(object):
    '''
    collected definitions to measure the length of an ellipse


    '''

    @staticmethod
    def ellip_drop(A,B,drop=0.4):
        '''
        given a list of axis lengths, calculate the length of the bar based on some specified ellipticity drop

        '''
        found = False
        j = 2
        while found==False:
            d = (1.-B[j-1]/A[j-1]) - (1.-B[j]/A[j])
            if d > drop:
                found = True
                print('INDEX VALUE is {0:d}'.format(j-1))
            j += 1
            if j==len(A):
                found = True
                j=2
        return A[j-2]



def ellip_drop_below(A,B,drop=0.4):
    '''
    where does the ellipticity first drop below some value?

    '''
    d = (1.-B/A)
    lessthan = np.where( d >= drop )[0]
    if len(lessthan) > 0:
        if np.max(lessthan) < 1:
            minbin = 1
        else:
            minbin = np.max(lessthan)
    else:
        minbin = 1
    return A[minbin]



def max_ellip_drop(A,B):
    edrop = np.ediff1d((1.-B/A),to_end=0.)
    return A[ np.where(np.min(edrop)==edrop)[0]]




#
# MUNOZ13 proposes several bar length metrics, reproduced here:
#
def max_ellip(A,B):
    e = (1.-B/A)
    return A[ np.where(np.max(e)==e)[0]]


def ellip_change(A,B,change=0.1):
    e = (1.-B/A)
    ellip_index = np.where(np.max(e)==e)[0]
    max_ellip_value = e[ellip_index]
    ellip_diff = 0.
    while (ellip_diff < change):
        ellip_index += 1
        ellip_diff = abs(e[ellip_index] - max_ellip_value)
    return A[ellip_index-1]


def pa_change(A,B,change=10.):
    # change must be in degrees
    e = (1.-B/A)
    pa = np.arctan(B/A)
    ellip_index = np.where(np.max(e)==e)[0]
    pa_value = pa[ellip_index]
    pa_diff = 0.
    while (pa_diff < change):
        ellip_index += 1
        pa_diff = abs(pa[ellip_index] - pa_value)*180./np.pi
    return A[ellip_index-1]
