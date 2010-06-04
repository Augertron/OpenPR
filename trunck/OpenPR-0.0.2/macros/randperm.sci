//generate a random permutation of the integers from 1 to n

function idx = randperm(n)
  
  number = rand(1, n);
  [tmp, idx] = gsort(number);
  
endfunction
