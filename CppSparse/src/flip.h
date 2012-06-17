//-*- mode: c++ -*-
#ifndef _FLIP_H_
#define _FLIP_H_
namespace
{
    template<class idx_type, bool sign>
    struct flip_t 
    {
    };

    template<class idx_type>
    struct flip_t<idx_type, true>
    {
        bool is_flipped(idx_type value)
        {
            return value < 0;
        }

        idx_type flip(idx_type value)
        {
            return -2 - value;
        }
    };

    template<class idx_type>
    struct flip_t<idx_type, false>
    {
        static const idx_type half = static_cast<idx_type>(-1)/2  + 1;

        static bool is_flipped(idx_type value)
        {
            return value > half;
        }

        static idx_type flip(idx_type value)
        {
            return value + half;
        }
    };

    template<class idx_type>
    inline bool is_flipped(idx_type val)
    {
        return flip_t<idx_type, std::is_signed<idx_type>::value>::is_flipped(val);
    }

    template<class idx_type>
    inline idx_type flip(idx_type val)
    {
        return flip_t<idx_type, std::is_signed<idx_type>::value>::flip(val);
    }

    template<class idx_type>
    inline idx_type unflip(idx_type val)
    {
        return is_flipped(val) ? flip(val) : val;
    }
}
#endif
