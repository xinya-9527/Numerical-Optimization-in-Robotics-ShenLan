#include "iostream"
#include <vector>
#include <math.h>

// define the dimension of optimization variables
#define  N   36
// calculate the vector norm
float norm(std::vector<float>& nums)
{
    float sum = 0.0f;
    for(std::vector<float>::iterator iter = nums.begin(); iter != nums.end(); iter++)
    {
        sum += *iter * *iter;
    }
    return sqrt(sum);
}
// calculate the gradient of rosenbrock
void cal_rosenbrock_grad(const std::vector<float>& nums, std::vector<float>& grad)
{
    for(uint8_t i = 0; i < nums.size()/2; i++)
    {
        grad[2*i] = 400 * (nums[2*i] * nums[2*i] - nums[2*i+1]) * nums[2*i] + 2 * (nums[2*i] - 1);
        grad[2*i+1] = -200 * (nums[2*i] * nums[2*i] - nums[2*i+1]);
    }
}
// calculate the value of rosenbrock
float rosenbrock(const std::vector<float>& nums)
{
    float sum = 0.0f;
    for(uint8_t i = 0; i < nums.size()/2; i++)
    {
        sum += 100 * (nums[2*i] * nums[2*i] - nums[2*i+1]) * (nums[2*i] * nums[2*i] - nums[2*i+1]);
        sum += (nums[2*i] - 1) * (nums[2*i] - 1);
    }
    return sum;
}
// operator overloading 
std::vector<float> operator*(const float& dist, const std::vector<float>& nums)
{
    std::vector<float> result(nums.size());
    for(uint8_t i = 0; i < nums.size(); i++)
    {
        result[i] = dist * nums[i];
    }
    return result;
}
std::vector<float> operator-(const std::vector<float>& numsa, const std::vector<float>& numsb)
{
    std::vector<float> result(numsa.size());
    for(uint8_t i = 0; i < numsa.size(); i++)
    {
        result[i] = numsa[i] - numsb[i];
    }
    return result;
}
std::ostream & operator<<(std::ostream& out, const std::vector<float>& nums)
{
    for(uint8_t i = 0; i < nums.size(); i++)
    {
        std::cout << nums[i] << "   ";
    }
    return out;
}

int main()
{
    std::vector<float> X(N);
    std::vector<float> gradient(N);
    if(X.size() % 2)
    {
        std::cout << "The dimension of X is incompatible." << std::endl;
        return 1;
    }
    // Initial state of optimization
    for(std::vector<float>::iterator iter = X.begin(); iter != X.end(); iter++)
    {
        *iter = float(rand()%1000)/100.0f;
        std::cout << *iter <<std::endl;
    }

    cal_rosenbrock_grad(X, gradient);
    float gradient_norm = norm(gradient);
    int cnt = 0;
    while (gradient_norm > 1e-3)
    {
        cnt ++;
        cal_rosenbrock_grad(X, gradient);
        float stepping = 1.0f;
        while (rosenbrock(X) < rosenbrock(X - stepping * gradient))
        {
            stepping *= 0.5f;
        }
        X = X - stepping * gradient;
        cal_rosenbrock_grad(X, gradient);
        gradient_norm = norm(gradient);
        std::cout << "iteration " << cnt << " times; gradient is " <<  gradient_norm << std::endl;
    }
    std::cout << "X: " << X << std::endl;
    return 0;
}