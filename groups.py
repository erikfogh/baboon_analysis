
import unittest

from gwf import Workflow, AnonymousTarget

gwf = Workflow()

#Changed it so it adds suffixes instead.

class Group:

    def __init__(self, workflow, suffix, sep="_"):
        self.workflow = workflow
        self.suffix = suffix
        self.sep = sep

    def _make_name(self, name):
        if name is None:
            return self.suffix
        return "{}{}{}".format(name, self.sep, self.suffix)

    def target(self, name, *args, **kwargs):
        suffixed_name = self._make_name(name)
        return self.workflow.target(*args, name=suffixed_name, **kwargs)

    def target_from_template(self, name, *args, **kwargs):
        suffixed_name = self._make_name(name)
        return self.workflow.target_from_template(suffixed_name, *args, **kwargs)

    def map(self, *args, name=None, **kwargs):
        if callable(name):
            def name_func(*args, **kwargs):
                return self._make_name(name(*args, **kwargs))
            return self.workflow.map(*args, name=name_func, **kwargs)
        else:
            return self.workflow.map(*args, name=self._make_name(name), **kwargs)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass


def my_template(x):
    return AnonymousTarget(inputs=[x], outputs=[], options={}, spec="")


class TestGroup(unittest.TestCase):
    def assertHasTargets(self, workflow, targets):
        self.assertEqual(set(workflow.targets.keys()), set(targets))

    def setUp(self):
        self.gwf = Workflow()

    def test_group_with_target_1(self):
        with Group(self.gwf, suffix="foo") as g:
            g.target("a", inputs=[], outputs=[]) << ""
            g.target("b", inputs=[], outputs=[]) << ""

        self.assertHasTargets(self.gwf, ["foo_a", "foo_b"])

    def test_group_with_target_2(self):
        with Group(self.gwf, suffix="bar") as g:
            g.target("a", inputs=[], outputs=[]) << ""
            g.target("b", inputs=[], outputs=[]) << ""

        self.assertHasTargets(self.gwf, ["bar_a", "bar_b"])

    def test_group_with_target_from_template(self):
        with Group(self.gwf, suffix="foobar") as g:
            g.target_from_template("baz", my_template("a"))

        self.assertHasTargets(self.gwf, ["foobar_baz"])

    def test_group_with_map_1(self):
        with Group(self.gwf, suffix="step_one") as g:
            lst = ["a", "b", "c", "d"]
            g.map(my_template, lst)

        self.assertHasTargets(self.gwf, ["step_one_0", "step_one_1", "step_one_2", "step_one_3"])

    def test_group_with_map_2(self):
        with Group(self.gwf, suffix="step_one") as g:
            lst = ["a", "b", "c", "d"]
            g.map(my_template, lst, name="baz")

        self.assertHasTargets(self.gwf, ["step_one_baz_0", "step_one_baz_1", "step_one_baz_2", "step_one_baz_3"])

    def test_group_with_map_nested(self):
        with Group(self.gwf, suffix="one") as g, Group(g, suffix="two") as g2:
            lst = ["a", "b"]
            g.map(my_template, lst, name="foo")

            lst = ["c", "d"]
            g2.map(my_template, lst, name="bar")

        self.assertHasTargets(self.gwf, ['one_foo_1', 'one_two_bar_0', 'one_foo_0', 'one_two_bar_1'])