#include "Piraha.hpp"

using namespace cctki_piraha;

bool NegLookAhead::match(std::shared_ptr<Matcher> m) {
    int pos = m->pos;
    bool b = pattern->match(m);
    m->pos = pos;
    return !b;
}
