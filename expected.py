


#idea is given any sequence of values, calculate the expected value 

def helper(p_array): #expected value sorted from high to low 
    expected = 0
    for i in range(len(p_array)-1, -1, -1):
        flipped = i + 1
        flipped *= p_array[i] if i != len(p_array)-1 else 1
        expected += flipped 
        expected *= (1- p_array[i - 1]) if i > 0 else 1

    #just get the heads version 
    return expected

     
def expected(i, arrayH, arrayT): #i is the probability of the first value 

    print(arrayH, arrayT)
    answer = 1
    answer += i * helper(arrayT) #use arrayT when you are trying to get tails
    answer += (1-i) * helper(arrayH) #use arrayH when you are trying to get heads
    return answer


array = [0.2, 0.9, 0.8, 0.7, 0.8] #array tracks the chance of getting heads 
chance_of_heads = 0.1

#flip it for the optimal value 
optimal = expected(chance_of_heads, sorted(array, reverse=True), [1-val for val in sorted(array, reverse=True)][::-1])
given = expected(chance_of_heads, array, [1-val for val in array])

print(optimal)
print(given)