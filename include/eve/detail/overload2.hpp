//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#define EVE_DISPATCH() delay_t{}, EVE_CURRENT_API{}

namespace eve::detail
{
  template<typename Tag, typename Caller> struct make_callable
  {
    using tag_type    = Tag;

    template<conditional_expr C, typename... Args>
    EVE_FORCEINLINE constexpr auto operator()(C&& c, Args &&... args) const noexcept
    {
      if constexpr( tag_dispatchable<tag_type,C,Args...> )
      {
        return tagged_dispatch(tag_type{}, std::forward<C>(c), std::forward<Args>(args)...);
      }
      else
      {
        return Caller::conditional_call(std::forward<C>(c), std::forward<Args>(args)...);
      }
    }

    template<typename... Args>
    EVE_FORCEINLINE constexpr auto operator()(Args &&... args) const noexcept
    {
      if constexpr( tag_dispatchable<tag_type,Args...> )
      {
        return tagged_dispatch(tag_type{}, std::forward<Args>(args)...);
      }
      else
      {
        return Caller::call(std::forward<Args>(args)...);
      }
    }
  };

  template<typename Tag, typename Caller> struct make_conditional
  {
    template<value Condition>
    EVE_FORCEINLINE constexpr auto operator[](Condition const &c) const noexcept
    {
      return  [cond = if_(to_logical(c))](auto const&... args) EVE_LAMBDA_FORCEINLINE
              {
                return make_callable<Tag,Caller>{}(cond, args...);
              };
    }

    template<conditional_expr Condition>
    EVE_FORCEINLINE constexpr auto operator[](Condition const &c) const noexcept
    {
      return  [c](auto const&... args) EVE_LAMBDA_FORCEINLINE
              {
                return make_callable<Tag,Caller>{}(c, args...);
              };
    }
  };
}

#define EVE_MAKE_CALLABLE2(TAG, NAME)                                                              \
inline std::ostream& operator<<(std::ostream& os, detail::callable_##TAG const&)                   \
{                                                                                                  \
  return os << #NAME;                                                                              \
}                                                                                                  \
inline detail::callable_##TAG const NAME = {}                                                      \
/**/

/*
namespace eve
{
  //================================================================================================
  // decorator definition, detection, combination and application to callables
  //================================================================================================
  struct decorator_ {};
  template<typename ID> concept decorator = std::derived_from<ID,decorator_>;

  template<typename Decoration> struct decorated;
  template<typename Decoration, typename... Args>
  struct decorated<Decoration(Args...)> : decorator_
  {
    using base_type = Decoration;

    template<decorator Decorator>
    constexpr EVE_FORCEINLINE auto operator()(Decorator d) const noexcept
    {
      return Decoration::combine(d);
    }

    template <typename Function>
    struct fwding_lamda
    {
      Function f;

      template <typename... X>
      constexpr EVE_FORCEINLINE auto operator()(X&&... x)
      {
        return f(decorated{}, std::forward<X>(x)...);
      }
    };

    template<typename Function>
    constexpr EVE_FORCEINLINE auto operator()(Function f) const noexcept
    {
      if constexpr( requires{ Decoration{}(f); } )  return Decoration{}(f);
      else                                          return fwding_lamda<Function>{f};
    }
  };

  namespace detail
  {
    template<typename C> struct if_;

    //==============================================================================================
    // basic type to support delayed calls
    struct delay_t {};

    //==============================================================================================
    // Extension point for centralizing asserts & static_asserts
    template<typename Tag, typename... Args>
    void check(delay_t const&, Tag const&, Args const&... ) {}

    //==============================================================================================
    // User-facing tag-dispatch helper
    template <typename Tag, typename... Args>
    concept tag_dispatchable = requires(Tag tag, Args&&... args)
    {
      { tagged_dispatch(tag, std::forward<Args>(args)...) };
    };
  }
}
*/
