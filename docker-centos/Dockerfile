FROM centos:7

MAINTAINER Stefan Seritan <sseritan@gmail.com>

# Repos
RUN yum install -y \
    centos-release-scl \
    https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
RUN yum-config-manager --enable rhel-server-rhscl-7-rpms
RUN curl -sL https://rpm.nodesource.com/setup_12.x | bash -

# Cmake3, GCC 7, Node.js 12
RUN yum install -y cmake3 devtoolset-7 nodejs rh-git218
RUN ln -s /usr/bin/cmake3 /usr/bin/cmake

# Manually enable devtoolset-7 and git 2.18
ENV PATH="/opt/rh/rh-git218/root/usr/bin:/opt/rh/devtoolset-7/root/usr/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/rh/httpd24/root/usr/lib64:/opt/rh/devtoolset-7/root/usr/lib64:/opt/rh/devtoolset-7/root/usr/lib:/opt/rh/devtoolset-7/root/usr/lib64/dyninst:/opt/rh/devtoolset-7/root/usr/lib/dyninst:/opt/rh/devtoolset-7/root/usr/lib64:/opt/rh/devtoolset-7/root/usr/lib"
