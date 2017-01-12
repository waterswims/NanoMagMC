#ifndef _SPIN
#define _SPIN

#include <vector>

using namespace std;

class spin_type
{
protected:
    bool zero_flag;
public:
    spin_type(){}
    ~spin_type(){}
    virtual void rand_spin(){}
    virtual void zero_spin(){}
    virtual void change_spin(){}
    virtual void change_spin(spin_type* in){}
    virtual int i_access(){return 0;}
    virtual vector<double>& spin_access(){}
    virtual bool is_zero(){return true;}
    virtual void print(){}
};

class ising_spin: public spin_type
{
private:
    int spin;
public:
    ising_spin();
    ising_spin(const ising_spin& other);
    ~ising_spin(){}
    void rand_spin();
    void zero_spin();
    void change_spin(){spin *= -1;}
    void change_spin(spin_type* in){spin *= -1;}
    int i_access();
    void set_spin(int in);
    bool is_zero();
    ising_spin& operator=(const ising_spin& other);
    template<class T> ising_spin& operator=(const T& other){exit(301);}
    void print();
};

class heis_spin: public spin_type
{
private:
    vector<double> spins;
public:
    heis_spin();
    heis_spin(const heis_spin& other);
    ~heis_spin(){}
    void rand_spin();
    void zero_spin();
    void change_spin(){this->rand_spin();}
    void change_spin(spin_type* in);
    vector<double>& spin_access(){return spins;}
    bool is_zero();
    heis_spin& operator=(const heis_spin& other);
    template<class T> heis_spin& operator=(const T& other){exit(301);}
    void print();
};

#endif
