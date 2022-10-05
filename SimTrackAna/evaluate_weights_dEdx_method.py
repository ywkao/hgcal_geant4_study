#!/usr/bin/env python2

dummy_weight = 0.0
weightsPerLayer_V16 = [dummy_weight, 5.55, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86, 58.63, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 60.7, 83.08, 83.08, 83.43, 83.61, 83.61, 83.61, 83.61, 83.61, 83.61, 83.61]
weightsPerLayer_V16_cee = [dummy_weight, 5.55, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 9.4, 12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86, 13.54, 12.86]

#--------------------------------------------------

def calcWeights(weightsPerLayer):
    #for wei in zip(weightsPerLayer[:], weightsPerLayer[1:] + [weightsPerLayer[-1]]): print wei
    res = [sum(wei)/2. for wei in zip(weightsPerLayer[:], weightsPerLayer[1:] + [weightsPerLayer[-1]])]; res[0] = dummy_weight; return res;   

def run_default():
    # MeV
    weights = calcWeights(weightsPerLayer_V16)
    print weightsPerLayer_V16
    print weights

def print_element(my_container):
    c = my_container

    # print data in rows
    if isinstance(c, list):
        output = ""
        for i, ele in enumerate(c):
            if i==0: output = "%5.2f" % ele
            else:    output += ", %5.2f" % ele
        print output

    if isinstance(c, dict):
        print "idx, dE/dx, set0, set1, set2"
        for i in range(27):
            #output = "%d, %5.2f, %5.2f, %5.2f, %5.2f" % (i, weightsPerLayer_V16[i], c["set0"][i], c["set1"][i], c["set2"][i])
            output = "%d" % i
            output += ", %5.2f" % weightsPerLayer_V16[i]
            output += ", %5.2f" % c["set0"][i]
            output += ", %5.2f" % c["set1"][i]
            output += ", %5.2f" % c["set2"][i]
            print output

#--------------------------------------------------

def create_weight(tag):
    # shorten CEE dE/dx weight
    w = weightsPerLayer_V16

    # prepare output
    output = []
    if tag=="set0":
        output = calcWeights(weightsPerLayer_V16[:29])

    if tag=="set1":
        for i in range(27):
            if i==0:
                output.append(dummy_weight)
            elif i==1:
                weight = (w[1] + w[2] + w[3]) / 2.
                output.append(weight)
            elif i%2==0:
                output.append(0.)
            elif i%2==1:
                weight = (w[i-1] + w[i] + w[i+1] + w[i+2]) / 2.
                output.append(weight)
            else:
                print "[ERROR] this is not expected in set1 configuration. please check."

    if tag=="set2":
        indices_layer_with_zero_weight = [2, 4, 6, 7, 9, 11, 13, 16, 18, 20, 22, 24, 26]
        indices_layer_with_expected_weight = [3, 10, 12, 17, 19, 21, 23, 25]
        for i in range(27):
            if i==0:
                output.append(dummy_weight)
            elif i==1:
                weight = (w[1] + w[2] + w[3]) / 2.
                output.append(weight)
            elif i==5:
                weight = (w[i-1] + w[i] + w[i+1] + w[i+2] + w[i+3]) / 2.
                output.append(weight)
            elif i==8:
                weight = (w[i-2] + w[i-1] + w[i] + w[i+1] + w[i+2]) / 2.
                output.append(weight)
            elif i==14:
                weight = (w[i-1] + w[i] + w[i+1]) / 2.
                output.append(weight)
            elif i==15:
                weight = (w[i] + w[i+1] + w[i+2]) / 2.
                output.append(weight)
            elif i in indices_layer_with_zero_weight:
                output.append(0.)
            elif i in indices_layer_with_expected_weight:
                weight = (w[i-1] + w[i] + w[i+1] + w[i+2]) / 2.
                output.append(weight)
            else:
                print "[ERROR] this is not expected in set1 configuration. please check."

    return output


def test():
    weights = {}
    weights["set0"] = create_weight("set0")
    weights["set1"] = create_weight("set1")
    weights["set2"] = create_weight("set2")

    print_element(weightsPerLayer_V16[:29])
    print_element(weights["set0"])
    print_element(weights["set1"])
    print_element(weights["set2"])
    print_element(weights)


if __name__ == "__main__":
    #run_default()
    test()
