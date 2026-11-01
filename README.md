# Estimating-Radius-from-Gamow-factor-in-alpha-decay
This program estimates the touching radius 'a' for Ra226 alpha decay to Rn222 from experimental lifetime value from NuDat3 using Newtons method. 

Used AI to generate formatting; analytical calculations were done by hand to simplify program. Initial radius' are assumed with Wood-Saxon model and the Gamow factor is calculated from WKB approximation for even-even nuclei. 

Results:

### Newton iterations

| n | x_n | g(x_n) |
|:-:|:---:|:------:|
| 0 | 1.5353423663E-01 | 2.0542563642E-02 |
| 1 | 1.6228311864E-01 | 3.3678686472E-04 |
| 2 | 1.6243135105E-01 | 9.1766095389E-08 |
| 3 | 1.6243139147E-01 | 6.7723604502E-15 |

### Results

| Quantity | Value |
|:---------|------:|
| x = a/b | 0.1624313915 |
| Residual f(x) − C | 1.1102E-16 |
| Touching radius a | 8.1536 fm |
| R_α (r₀ = 1.2 fm) | 1.9049 fm |
| R(²²²Rn) = a − R_α | 6.2487 fm |
| Fitted r₀ from a | 11.1915 fm |
| Compare 1.2 × 222^(1/3) | 7.2661 fm |

## Future considerations

Extending the program to use FTCS-method and for even-odd and odd-odd nuclei could be an interesting future prospect to implement. 
