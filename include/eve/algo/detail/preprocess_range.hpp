//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/algo/concepts/eve_iterator.hpp>
#include <eve/algo/detail/convert.hpp>
#include <eve/algo/traits.hpp>

#include <iterator>
#include <type_traits>
#include <utility>

namespace eve::algo
{
  template <typename T, typename Cardinal>
  struct aligned_ptr_iterator;

  template <typename T, typename Cardinal>
  struct unaligned_ptr_iterator;

  namespace detail
  {
    template <typename Traits, typename I, typename S, typename ToOutput>
    struct preprocess_range_result
    {
     private:
      Traits traits_;
      I f_;
      S l_;
      ToOutput to_output_;
     public:
      preprocess_range_result(Traits traits, I f, S l, ToOutput to_output):
        traits_(traits),
        f_(f),
        l_(l),
        to_output_(to_output)
      {}

      Traits traits() const { return traits_; }

      I begin() const { return f_; }
      S end()   const { return l_; }

      template <typename I_>
      EVE_FORCEINLINE auto to_output_iterator(I_ it) const { return to_output_(it); }
    };

    template <typename Traits, typename I, typename S, typename ToOutput, typename Update>
    EVE_FORCEINLINE auto enhance_to_output(preprocess_range_result<Traits, I, S, ToOutput> prev, Update update )
    {
      return preprocess_range_result(
        prev.traits(),
        prev.begin(),
        prev.end(),
        [update, prev](auto it) { return update(prev.to_output_iterator(it)); }
      );
    }

    // Base case. Should validate that I, S are a valid iterator pair
    template <typename Traits, iterator I, sentinel_for<I> S>
    EVE_FORCEINLINE auto preprocess_eve_it_sentinel(Traits traits_, I f, S l)
    {
      if constexpr ( !std::same_as<typename I::value_type, iteration_type_t<Traits, I>> )
      {
        using T = iteration_type_t<Traits, I>;
        auto f_ = convert(f, eve::as<T>{});
        auto l_ = convert(l, eve::as<T>{});
        return detail::enhance_to_output(
          preprocess_eve_it_sentinel(traits_, f_, l_),
          [](unaligned_t<decltype(f_)> i)
          {
            return convert(i, eve::as<typename I::value_type>{});
          }
        );
      }
      else if constexpr ( typename I::cardinal{}() > forced_cardinal_t<Traits, typename I::value_type>{}() )
      {
        using N = forced_cardinal_t<Traits, typename I::value_type>;
        auto f_ = f.cardinal_cast(N{});
        return detail::enhance_to_output(
          preprocess_eve_it_sentinel(traits_, f_, l.cardinal_cast(N{})),
          [](unaligned_t<decltype(f_)> i)
          {
            return i.cardinal_cast(typename I::cardinal{});
          });
      }
      else
      {
        auto deduced = [] {
          if constexpr (partially_aligned_iterator<I>)
          {
            if constexpr ( std::same_as<I, S> && !always_aligned_iterator<I> ) return algo::traits(no_aligning, divisible_by_cardinal);
            else                                                               return algo::traits(no_aligning);
          }
          else
          {
            return algo::traits();
          }
        }();

        return detail::preprocess_range_result{
          default_to(traits_, deduced), f, l,
          [](unaligned_t<I> i) { return i; }
        };
      }
    }
  }

  struct preprocess_range_
  {
    template <typename Traits, std::contiguous_iterator I, typename S>
    EVE_FORCEINLINE auto operator()(Traits traits_, I f, S l) const
    {
      using T = std::remove_reference_t<decltype(*f)>;
      using it = unaligned_ptr_iterator<T, forced_cardinal_t<Traits, typename std::iterator_traits<I>::value_type>>;

      T* raw_f = nullptr;
      T* raw_l = raw_f;

      if (f != l)
      {
        raw_f = &*f;
        raw_l = raw_f + (l - f);
      }

      return detail::enhance_to_output(
        operator()(traits_, it{raw_f}, it{raw_l}),
        [f, raw_f](it i) { return f + (i.ptr - raw_f); }
      );
    }

    template <typename Traits, typename I, typename S>
      requires eve::detail::tag_dispatchable<preprocess_range_, Traits, I, S>
    EVE_FORCEINLINE auto operator()(Traits traits, I f, S l) const {
      return tagged_dispatch(*this, traits, f, l);
    }

    template <typename Traits, typename Rng>
      requires eve::detail::tag_dispatchable<preprocess_range_, Traits, Rng>
    EVE_FORCEINLINE auto operator()(Traits traits_, Rng&& rng) const {
      return tagged_dispatch(*this, traits_, std::forward<Rng>(rng));
    }

    template <typename Traits, typename Rng>
    EVE_FORCEINLINE auto operator()(Traits traits_, Rng&& rng) const {
      return operator()(traits_, rng.begin(), rng.end());
    }

    template <typename Traits, typename T, typename A>
    EVE_FORCEINLINE auto operator()(Traits traits_, eve::aligned_ptr<T, A> f, T* l) const
    {
      using N = forced_cardinal_t<Traits, T>;

      if constexpr (N{}() > A{}()) return operator()(traits_, f.get(), l);
      else
      {
        using aligned_it   = aligned_ptr_iterator<T, N>;
        using unaligned_it = unaligned_ptr_iterator<T, N>;

        return detail::enhance_to_output(
          operator()(traits_, aligned_it(f), unaligned_it(l)),
          [](unaligned_it i) { return i.ptr; }
        );
      }
    }

    template <typename Traits, typename T, typename A1, typename A2>
    EVE_FORCEINLINE auto operator()(Traits traits_, eve::aligned_ptr<T, A1> f, eve::aligned_ptr<T, A2> l) const
    {
      using N = forced_cardinal_t<Traits, T>;

           if constexpr ( N{}() > A2{}()) return operator()(traits_, f, l.get());
      else if constexpr ( N{}() > A1{}()) return operator()(traits_, f.get(), l);
      else
      {
        using aligned_it   = aligned_ptr_iterator<T, forced_cardinal_t<Traits, T>>;

        return detail::enhance_to_output(
          operator()(traits_, aligned_it(f), aligned_it(l)),
          [](unaligned_t<aligned_it> i) { return i.ptr; }
        );
      }
    }

    // Base case. Should validate that I, S are a valid iterator pair
    template <typename Traits, iterator I, sentinel_for<I> S>
    EVE_FORCEINLINE auto operator()(Traits traits_, I f, S l) const
    {
      return detail::preprocess_eve_it_sentinel(traits_, f, l);
    }

  } inline constexpr preprocess_range;
}