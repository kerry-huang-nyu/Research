


#idea is given any sequence of values, calculate the expected value 
def helper(p_array): #expect the to be sorted from low to high 
    expected = 0
    for i in range(len(p_array) , -1, -1):
        flipped = i 
        flipped *= p_array[i] if i != len(p_array) else 1
        expected += flipped 
        expected *= (1- p_array[i - 1]) if i != 0 else 1

    #just get the heads version 
    return expected

     
def expected(i, arrayH, arrayT): #i is the probability of the first value 
    answer = 0
    answer += i * helper(arrayH) #chance of getting heads 
    answer += (1-i) * helper(arrayT) #chance of getting tails 
    return answer


array = [0.2, 0.9, 0.8, 0.7, 0.8]


optimal = expected(0.1, sorted(array), [1-val for val in sorted(array)][::1])
given = expected(0.1, array, [1-val for val in array][::1])

print(optimal)
print(given)