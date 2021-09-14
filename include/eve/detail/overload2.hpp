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
    friend std::ostream& operator<<(std::ostream& os, make_callable const&)
    {
      return os << Caller::name;
    }

    using tag_type    = Tag;

    template<typename A0, typename... Args>
    EVE_FORCEINLINE constexpr auto operator()(A0&& a0, Args &&... args) const noexcept
    {
      if constexpr( tag_dispatchable<tag_type,A0,Args...> )
      {
        return tagged_dispatch(tag_type{}, std::forward<A0>(a0), std::forward<Args>(args)...);
      }
      else if constexpr( conditional_expr<std::remove_cvref_t<A0>> )
      {
        return Caller::conditional_call(std::forward<A0>(a0), std::forward<Args>(args)...);
      }
      else
      {
        return Caller::call(std::forward<A0>(a0), std::forward<Args>(args)...);
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

namespace eve::detail
{
  //================================================================================================
  // decorator definition & application to callables
  //================================================================================================
  template<typename Caller, typename Decorator> struct decorate
  {
    using type = Caller;
  };

/*
  struct decorator_ {};
  template<typename ID> concept decorator = std::derived_from<ID,decorator_>;
*/

  template<typename C> struct if_;
/*
    //==============================================================================================
    // basic type to support delayed calls
    struct delay_t {};

    //==============================================================================================
    // User-facing tag-dispatch helper
    template <typename Tag, typename... Args>
    concept tag_dispatchable = requires(Tag tag, Args&&... args)
    {
      { tagged_dispatch(tag, std::forward<Args>(args)...) };
    };
*/
}
