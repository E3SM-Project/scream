#include "expression_parser.hpp"

#include <Kokkos_Core.hpp>

namespace scream {

// enum class BinOp {
//   Add,
//   Sub,
//   Prod,
//   Div,
//   Exp,
// };

// enum class UnaryOp {
//   Minus,
//   Parens
// };

// enum class Op {
//   Add,
//   Sub,
//   Mult,
//   Div,
//   Exp,
//   Minus,
//   Load
// };

// template<Op op1, Op op2>
// struct Precedes;

// // template<typename S1, typename S2>
// // struct Promote {
// // };

// // template<typename S,int N>
// // struct Promote<S,ekat::Pack<S,N>> {
// //   using type = ekat::Pack<S,N>;
// // };

// // template<typename S>
// // struct Promote<S,S> {
// //   using type = S;
// // };

// template<typename ScalarT,BinOp op>
// struct OpFunctor
// {
//   KOKKOS_FORCEINLINE_FUNCTION
//   static ScalarT apply (const ScalarT lhs, const ScalarT rhs);
// };

// template<typename ScalarT>
// struct OpFunctor<ScalarT,BinOp::Add>
// {
//   KOKKOS_FORCEINLINE_FUNCTION
//   static ScalarT apply (const ScalarT lhs, const ScalarT rhs)
//   {
//     return lhs+rhs;
//   }
// };
// template<typename ScalarT>
// struct OpFunctor<ScalarT,BinOp::Sub>
// {
//   KOKKOS_FORCEINLINE_FUNCTION
//   static ScalarT apply (const ScalarT lhs, const ScalarT rhs)
//   {
//     return lhs-rhs;
//   }
// };

// template<typename ScalarT, typename Node>
// struct EvalNode {
//   ScalarT lhs;
//   ScalarT rhs;
//   Op op;
//   ScalarT eval;

//   KOKKOS_FORCEINLINE_FUNCTION
//   void apply () {
//     switch (op) {
//       case Op::Add:   eval=lhs+rhs;       break;
//       case Op::Sub:   eval=lhs-rhs;       break;
//       case Op::Mult:  eval=lhs*rhs;       break;
//       case Op::Div:   eval=lhs/rhs;       break;
//       case Op::Exp:   eval=pow(lhs,rhs);  break;
//       case Op::Minus: eval=-lhs;          break;
//       case Op::Load:  eval=lhs;           break;
//       default:
//         EKAT_KERNEL_ERROR_MSG ("Unsupported/unrecognized operation.\n");
//     }
//   }
// };

// // Shortcut
// template<typename... ScalarT>
// using same_t = typename std::enable_if
//   <
//     ekat::SameType<ScalarT...>::value,
//     typename ekat::SameType<ScalarT...>::type
//   >::type;

// template<typename Op1, typename Op2>
// struct Compose {
//   template<typename... ScalarT>
//   KOKKOS_FORCEINLINE_FUNCTION
//   static same_t<ScalarT...> apply (const ScalarT... args);

//   template<typename ScalarT>
//   KOKKOS_FORCEINLINE_FUNCTION
//   ScalarT apply(const ScalarT arg0, const ScalarT arg2, const ScalarT arg3);
// };

ExprParser::
ExprParser  (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic (comm,params)
{
  m_expr = m_params.get<std::string>("expression");
  ekat::strip(m_expr,' ');
  validate_expression();
  // auto tree = ekat::parse_nested_list_no_sep(m_expr);

  // std::vector<EvalNode<Real>> nodes;
  // auto add_eval_nodes = [](const ekat::ParameterList& pl,auto&& recurse) {
  //   const auto ops_chars = std::vector<char>{'+','-','*','/','^'};
  //   int n = pl.get<int>("Num Entries");
  //   for (int i=0; i<n; ++i) {
  //     auto is_list = pl.get<std::string>(ekat::strint("Type",i))=="List";
  //     auto ei = ekat::strint("Entry",i);
  //     if (is_list) {
  //       recurse (pl.sublist(ei));
  //     } else {
  //       const auto& s = pl.get<std::string>(ei);
  //       auto tokens_seps = ekat::split(s,ops_chars);
  //       auto& tokens = tokens_seps.first;
  //       auto& seps = tokens_seps.second;
  //       EKAT_REQUIRE_MSG (i>0 || tokens.front().size()>0,"missing leading LHS\n");
  //       EKAT_REQUIRE_MSG (i<(n-1) || tokens.back().size()>0,"missing trailing RHS\n");
  //       int nt = tokens.size();
  //       // Sanity check
  //       for (int it=0; it<nt; ++it) {
  //         EKAT_REQUIRE_MSG (tokens[it].size()>0 || it==0 || it==(nt-1),
  //             "empty middle token\n");
  //       }

  //       while (nt>)
  //     }
  //   }
  // };
}

void ExprParser::compute_diagnostic_impl ()
{
  // int n = m_diagnostic_output.get_header().get_identifier().get_layout().size();
  // Kokkos::RangePolicy<DefaultDevice::execution_space> policy (0,n);

  // Kokkos::parallel_for (m_expr,policy,
  //                       KOKKOS_LAMBDA(const int i) {

  // });

}

void ExprParser::validate_expression () const
{
  std::vector<char> delims = {'*', '+', '/', '-', '^', '(', ')'};
  int open = 0;
  auto pos = m_expr.find_first_of(delims.data());
  while (pos!=std::string::npos) {
    auto c = m_expr[pos];
    if (c==')') {
      EKAT_REQUIRE_MSG(open>0,
          "[ExprParser] Error! Invalid expression.\n"
          "  - expression: " + m_expr + "\n"
          "  - reason    : invalid open/close parentheses sequence\n");
      --open;
    } else if (c=='(') {
      EKAT_REQUIRE_MSG (pos==0 or m_expr[pos-1]!=')',
          "[ExprParser] Error! Invalid expression.\n"
          "  - expression: " + m_expr + "\n"
          "  - reason    : opening parenthesis cannot follow closing parenthesis\n");
    } else {
      EKAT_REQUIRE_MSG (pos!=(m_expr.size()-1),
          "[ExprParser] Error! Invalid expression.\n"
          "  - expression: " + m_expr + "\n"
          "  - reason    : expression ends with arithmetic op character\n");
      if (pos==0) {
        EKAT_REQUIRE_MSG(c=='-',
            "[ExprParser] Error! Invalid expression.\n"
            "  - expression: " + m_expr + "\n"
            "  - reason    : expression cannot begin with '*', '+', or '^'\n");
      } else {
        auto p = m_expr[pos-1];
        EKAT_REQUIRE_MSG(p!='*' and p!='+' and p!='/' and p!='-' and p!='^',
            "[ExprParser] Error! Invalid expression.\n"
            "  - expression: " + m_expr + "\n"
            "  - reason    : cannot have two arithmetic operators in a row\n");

      }
    }
    pos = m_expr.find_first_of(delims.data(),pos+1);
  }
  EKAT_REQUIRE_MSG (open==0,
      "[ExprParser] Error! Invalid expression.\n"
      "  - expression: " + m_expr + "\n"
      "  - reason    : unmatched opening parenthesis\n");
}

} // namespace scream
