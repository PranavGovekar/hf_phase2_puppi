#ifndef HF_TOWERS_H
#define HF_TOWERS_H

#include "ap_int.h"

class sptower
{
public:
    ap_uint<14> energy;
    ap_uint<2> fb;

    sptower()
    {
        energy = 0;
        fb = 0;
    }

    ap_uint<10> gettower(void)
    {
        ap_uint<10> data;
        data  =
            ((ap_uint<10>)energy & 0xFF) |
            ((ap_uint<10>)fb<<8 & 0x300) ;
        return data ;
    }

    void fillhftower(ap_uint<10> i)
    {
        this->energy = i.range(7, 0);
        this->fb = i.range(9, 8);
    }

    sptower(ap_uint<10> i)
    {
        this->energy = i.range(7, 0);
        this->fb = i.range(9, 8);
    }

    sptower(const sptower& rhs)
    {
        energy=rhs.energy;
        fb=rhs.fb;
    }

    sptower& operator=(const sptower& rhs)
    {
        this->energy=rhs.energy;
        this->fb=rhs.fb;
        return *this;
    }

    ap_uint<8> Energy(void)
    {
        return energy;
    }
    ap_uint<2> Fb(void)
    {
        return fb;
    }
};


class hftower {
public:
    ap_uint<10> energy;
    ap_uint<2> fb;
    ap_uint<8> phi;
    ap_uint<5> eta;

    hftower() : energy(0), fb(0), eta(0), phi(0) {}

    hftower(ap_uint<10> i) {
        this->energy = i.range(7, 0);
        this->fb = i.range(9, 8);
        eta = 0;
        phi = 0;
    }

    hftower(const hftower& rhs) {
        energy = rhs.energy;
        fb = rhs.fb;
        eta = rhs.eta;
        phi = rhs.phi;
    }

    hftower& operator=(const hftower& rhs) {
        energy = rhs.energy;
        fb = rhs.fb;
        eta = rhs.eta;
        phi = rhs.phi;
        return *this;
    }

    ap_uint<10> gettower(void) const {
        ap_uint<10> data =
            ((ap_uint<10>)energy & 0xFF) |
            ((ap_uint<10>)fb << 8 & 0x300);
        return data;
    }

    void fillhftower(ap_uint<10> i) {
        energy = i.range(7, 0);
        fb = i.range(9, 8);
    }

    ap_uint<8> Energy(void) const { return energy; }
    ap_uint<2> Fb(void) const { return fb; }
};


#endif
