function [element,product,summation]=matrixops(mat1,mat2)
% [element,product,summation]=MATRIXOPS(mat1,mat2) is a function which takes 
% as inputs two marices (must be of the same size) and outputs the 
% element-by-element product, the product mat1'*mat2, and the summation
% of all elements of mat1 and mat2
if size(mat1)~=size(mat2)
    error('mat1 and mat2 must have the same size')
end
element = mat1.*mat2;
product = mat1'*mat2;
summation = sum(sum(mat1+mat2));