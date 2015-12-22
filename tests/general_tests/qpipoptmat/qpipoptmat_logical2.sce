//Find x in R^6 such that:
// A simple example with constraints

Aeq= [1,-1,1,0,3,1;
-1,0,-3,-4,5,6;
2,5,3,0,1,0];
beq=[1; 2; 3];
A= [0,1,0,1,2,-1;
-1,0,2,1,1,0];
b = [-1; 2.5];
lb=[-1000; -10000; 0; -1000; -1000; -1000];
ub=[10000; 100; 1.5; 100; 100; 1000];
param = list("MaxIter", 300, "CpuTime",100);
//and minimize 0.5*x'*Q*x + p'*x with
f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
[xopt,fopt,exitflag,output,lambda]=qpipoptmat(H,f,A,b,Aeq,beq,lb,ub,[],param)
disp("xopt",xopt,"fopt",fopt,"exitflag",exitflag,"output",output,"lambda",lambda)

//Output
//
//Optimal Solution Found.
//
//   lower: [1x6 constant]
//   upper: [1x6 constant]
//   ineqlin: [0.3091368,1.197D-12]
//   eqlin: [-1.5564027,-0.1698164,-0.7054782]
// 
// lambda   
// 
//   Iterations: 13
// 
// output   
// 
//  0  
// 
// exitflag   
// 
//  - 14.843248  
// 
// fopt   
// 
//    1.7975426  
//  - 0.3381487  
//    0.1633880  
//  - 4.9884023  
//    0.6054943  
//  - 3.1155623  
// 
// xopt   
