function gaussian_kernel(a, b) result(kernel)

    double precision, dimension(:), intent(in) :: a
    double precision, dimension(:), intent(in) :: b

    double precision :: kernel
    integer :: i
    double precision :: temp
    double precision :: sigma
   
    kernel = 0.0d0 
    sigma = 724.0d0

    do i = 1, size(a)
        temp = a(i) - b(i)
        kernel = kernel + temp * temp
    enddo

    kernel = exp(-0.5d0 * sqrt(kernel) / (sigma*sigma))

end function gaussian_kernel


function laplace_kernel(a, b) result(kernel)

    double precision, dimension(:), intent(in) :: a
    double precision, dimension(:), intent(in) :: b

    double precision :: kernel

    kernel = sum(abs(a(:) - b(:)))
    kernel = exp(-0.00025118864315095795d0 * kernel)

end function laplace_kernel

subroutine kgaussian_kernel(a, na, b, nb, k, sigma)

    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b

    integer, intent(in) :: na, nb

    double precision, dimension(:,:), intent(inout) :: k
    double precision, intent(in) :: sigma

    double precision, allocatable, dimension(:) :: temp

    double precision :: inv_sigma
    integer :: l

    inv_sigma = -0.5d0 / (sigma*sigma)

    allocate(temp(size(a, dim=1)))

!$OMP PARALLEL DO PRIVATE(temp)
    do i = 1, nb 
        do j = 1, na
            temp(:) = a(:,j) - b(:,i)
            K(j,i) = exp(inv_sigma * sqrt(sum(temp*temp)))
        enddo
    enddo
!$OMP END PARALLEL DO

    deallocate(temp)
end subroutine kgaussian_kernel

subroutine klaplace_kernel(a, na, b, nb, k, sigma)

    double precision, dimension(:,:), intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b

    integer, intent(in) :: na, nb

    double precision, dimension(:,:), intent(inout) :: k
    double precision, intent(in) :: sigma

    double precision :: inv_sigma

    inv_sigma = -1.0d0 / sigma

!$OMP PARALLEL DO
    do i = 1, nb 
        do j = 1, na
        k(j,i) = exp(inv_sigma * sum(abs(a(:,j) - b(:,i))))
        enddo
    enddo
!$OMP END PARALLEL DO

end subroutine klaplace_kernel
