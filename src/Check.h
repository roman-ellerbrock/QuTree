//
// Created by Roman Ellerbrock on 11/18/21.
//

#ifndef CHECK_H
#define CHECK_H

#define CHECK(X) ({int __val = (X); __val == -1 ? \
({ fprintf(stderr, "ERROR (" __FILE__ ":%d) -- %s\n",__LINE__,strerror(errno)); \
exit(-1);-1;}) : __val; })


#endif //CHECK_H
