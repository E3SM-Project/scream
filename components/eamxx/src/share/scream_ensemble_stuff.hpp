/*
f_layout: (ncols,ndims,nlevs)

f_storage_layout: (ncols,ndims,nlevs*ens_size)

template<typename R,typename T>
PackVectorReduce {
  template<typename Pack>
  void reduce (const T& t, const Pack& p) {
    R r(p);
    t.vector_reduce(r);
  }
}

template<typename R>
PackVectorReduce<R,HostThreadTeamMember>
{
  template<typename Pack>
  void reduce (const T&, const Pack&) {}
};

template<typename ExecSpace>
Kokkos::TeamPolicy<ExecSpace>::member_type
get_empty_team_member ();

template<>
Kokkos::TeamPolicy<Kokkos::Cuda>::member_type
get_empty_team_member () {
  return Kokkos::Impl::CudaTeamMember(nullptr,0,0,nullptr,0,0,0);
}



op+(Pack p1, Pack p2) {
  auto team = get_empty_team_member();
  auto rang = Kokkos::ThreadVectorRange(team,Pack::n);
  Pack r;
  Kokkos::parallel_for(range,[&](int i){
    r[i] = p1[i]+p2[i];
  });
  vector_reduce<Kokkos::Sum,TeamMember>(r);
}

*/
