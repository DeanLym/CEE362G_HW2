function [] = hw2()
    data = load('2015Assign1_1.txt');
    t = data(:,1);
    y = data(:,2);
    % number of unknowns and number of observations
    m = 500; n = 150;
    n = 200;
    H = shaw(n); 
    X = ones(n,1);
    q1();
end


function [] = q1()





end


