hw06_worker = hw06();

hw_assert(abs ( hw06_worker.p1(@f1_test, 1, 1) - 1/4) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f1_test, 4, 1) - 0.328125) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f2_test, 1, 1) - 1/8) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f2_test, 4, 1) - 0.2421875) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f3_test, 1, 1) - 1/16) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f3_test, 4, 1) - 0.189697265625) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f4_test, 1, 1) - 1/32) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f4_test, 4, 1) - 0.15393066406250) < 1e-8);

hw_assert(abs ( hw06_worker.p1(@f1_test, 1, 2) - 1/2) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f1_test, 4, 2) - 0.34375) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f2_test, 1, 2) - 1/2) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f2_test, 4, 2) - 0.265625) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f3_test, 1, 2) - 1/2) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f3_test, 4, 2) -  0.220703125) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f4_test, 1, 2) - 1/2) < 1e-8);
hw_assert(abs ( hw06_worker.p1(@f4_test, 4, 2) - 0.1923828125) < 1e-8);

hw_assert(abs ( hw06_worker.p1(@f1_test, 10, 3) - 1/3) < 1e-5);
hw_assert(abs ( hw06_worker.p1(@f2_test, 10, 3) - 1/4) < 1e-5);
hw_assert(abs ( hw06_worker.p1(@f3_test, 10, 3) - 1/5) < 1e-5);
hw_assert(abs ( hw06_worker.p1(@f4_test, 10, 3) - 1/6) < 1e-5);
hw_assert(abs ( hw06_worker.p1(@f5_test, 10, 3) - 1/7) < 1e-5);
hw_assert(abs ( hw06_worker.p1(@f6_test, 10, 3) - (1-cos(1))) < 1e-5);
hw_assert(abs ( hw06_worker.p1(@f7_test, 10, 3) - (exp(1)-1)) < 1e-5);
hw_assert(abs ( hw06_worker.p1(@f8_test, 10, 3) - pi/4) < 1e-5);

function ret = f1_test(x)
    ret = x.^2;
end

function ret = f2_test(x)
    ret = x.^3;
end

function ret = f3_test(x)
    ret = x.^4;
end

function ret = f4_test(x)
    ret = x.^5;
end

function ret = f5_test(x)
    ret = x.^6;
end

function ret = f6_test(x)
    ret = sin(x);
end

function ret = f7_test(x)
    ret = exp(x);
end

function ret = f8_test(x)
    ret = 1./(1 + x.^2);
end

function hw_assert(X)
    if X; fprintf('\t PASS\n'); else; fprintf('\t FAIL\n'); end
end