#pragma once

#include <type_traits>
#include <utility>

namespace utl {

namespace {
  
// type trait to get the type of the first argument of a function
template <typename Return, typename Arg0, typename... Args>
Arg0 first_argument_impl(Return (*)(Arg0, Args...));

template <typename Return, typename Function, typename Arg0, typename... Args>
Arg0 first_argument_impl(Return (Function::*)(Arg0, Args...));

template <typename Return, typename Function, typename Arg0, typename... Args>
Arg0 first_argument_impl(Return (Function::*)(Arg0, Args...) const);

template <typename Function>
decltype(first_argument_impl(&Function::operator())) first_argument_impl(Function);

}
// when passing T to the template alias the type will be the return type of first_argument_impl called with a T object
template <typename T>
using first_argument = decltype(first_argument_impl(std::declval<T>()));

}  // namespace utl